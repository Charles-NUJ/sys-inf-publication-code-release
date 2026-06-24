include(joinpath(@__DIR__, "interference_all_to_all_demo.jl"))
include(joinpath(@__DIR__, "systematic_dfs_baseline.jl"))

struct ComparisonResult
    method::String
    iteration::Int
    period::Int
    delay_sum::Int
    average_aoi::Float64
    runtime::Float64
    mean_parallelism::Float64
    maximum_parallelism::Int
    schedule::Vector{Vector{ScheduledBroadcast}}
    generation_times::Vector{Int}
    arrival_times::Matrix{Int}
    first_parent::Matrix{Int}
    first_reception_slot::Matrix{Int}
end

function transmitter_usage(schedule, node_count)
    usage = zeros(Float64, node_count)
    for actions in schedule, action in actions
        usage[action.transmitter] += 1
    end
    return usage
end

function run_one_schedule(graph, conflict_graph, method, iteration,
                          rng, noise_scale, node_prices;
                          score_parameters=DEFAULT_SCORE_PARAMETERS,
                          maximum_transmitters=typemax(Int))
    policy = method == "Coverage-first" ? :coverage :
             method == "Age-first" ? :age : :dual
    output = nothing
    runtime = @elapsed output = schedule_all_to_all(
        graph, conflict_graph; policy=policy, rng=rng,
        noise_scale=noise_scale, node_prices=node_prices,
        score_parameters=score_parameters,
        maximum_transmitters=maximum_transmitters)
    schedule, generation_times, arrival_times,
    first_parent, first_reception_slot = output
    validate_schedule(schedule, conflict_graph)
    period, delay_sum, average_aoi = evaluate_schedule(
        schedule, generation_times, arrival_times)
    parallelism = [length(actions) for actions in schedule]
    return ComparisonResult(
        method, iteration, period, delay_sum, average_aoi, runtime,
        sum(parallelism) / period, maximum(parallelism), schedule,
        generation_times, arrival_times, first_parent, first_reception_slot)
end

"""
Run repeated randomized construction for one policy.

For the dual-guided policy, prices are updated from node transmission loads.
An above-average node load raises its price, encouraging subsequent iterations
to use alternative relays and less congested spatial regions.
"""
function run_method(graph, conflict_graph, method, iterations, seed;
                    perturbation_scale=0.65,
                    dual_perturbation_scale=0.50,
                    dual_step_scale=0.08,
                    score_parameters=DEFAULT_SCORE_PARAMETERS,
                    maximum_transmitters=typemax(Int))
    node_count = nv(graph)
    prices = zeros(Float64, node_count)
    incumbent = nothing
    all_results = ComparisonResult[]
    incumbent_curve = Float64[]
    current_curve = Float64[]
    method_offset = method == "Coverage-first" ? 101 :
                    method == "Age-first" ? 202 : 303

    for iteration in 1:iterations
        rng = MersenneTwister(seed + 10_000 * iteration + method_offset)
        noise_scale = iteration == 1 ? 0.0 :
                      method == "Dual-guided" ? dual_perturbation_scale :
                      perturbation_scale
        result = run_one_schedule(graph, conflict_graph, method, iteration,
                                  rng, noise_scale, prices;
                                  score_parameters=score_parameters,
                                  maximum_transmitters=maximum_transmitters)
        push!(all_results, result)
        push!(current_curve, result.average_aoi)
        if incumbent === nothing || result.average_aoi < incumbent.average_aoi
            incumbent = result
        end
        push!(incumbent_curve, incumbent.average_aoi)

        if method == "Dual-guided"
            usage = transmitter_usage(result.schedule, node_count)
            centered_load = usage .- sum(usage) / node_count
            step_size = dual_step_scale / sqrt(iteration)
            prices .= max.(0.0, prices .+ step_size .* centered_load)
        end
    end
    return incumbent, all_results, incumbent_curve, current_curve
end

function save_comparison_csv(path, results)
    open(path, "w") do output
        println(output, "method,iteration,period,delay_sum,average_aoi,runtime," *
                        "mean_parallelism,maximum_parallelism")
        for result in results
            println(output, "$(result.method),$(result.iteration),$(result.period)," *
                            "$(result.delay_sum),$(result.average_aoi),$(result.runtime)," *
                            "$(result.mean_parallelism),$(result.maximum_parallelism)")
        end
    end
end

function complement_dsatur_coloring_upper_bound(conflict_graph)
    node_count = nv(conflict_graph)
    colors = zeros(Int, node_count)
    complement_degrees = [node_count - 1 - degree(conflict_graph, node)
                          for node in 1:node_count]
    for _ in 1:node_count
        uncolored = findall(==(0), colors)
        selected = uncolored[1]
        best_key = (-1, -1)
        for node in uncolored
            neighbor_colors = Set{Int}()
            for other in 1:node_count
                if other != node && !has_edge(conflict_graph, node, other) &&
                   colors[other] != 0
                    push!(neighbor_colors, colors[other])
                end
            end
            key = (length(neighbor_colors), complement_degrees[node])
            if key > best_key
                selected, best_key = node, key
            end
        end
        forbidden = Set(colors[other] for other in 1:node_count
                        if other != selected &&
                           !has_edge(conflict_graph, selected, other) &&
                           colors[other] != 0)
        color = 1
        while color in forbidden
            color += 1
        end
        colors[selected] = color
    end
    # alpha(G) = omega(complement(G)) <= chi(complement(G)) <= this coloring.
    return maximum(colors)
end

function common_aoi_lower_bound(graph, conflict_graph)
    node_count = nv(graph)
    maximum_degree = maximum(degree(graph))
    per_source_broadcast_lb = cld(node_count - 1, maximum_degree)
    total_broadcast_lb = node_count * per_source_broadcast_lb
    independence_upper_bound = complement_dsatur_coloring_upper_bound(conflict_graph)
    period_lb = ceil(Int, total_broadcast_lb / independence_upper_bound - 1e-9)
    delay_lb = distance_delay_lower_bound(graph)
    aoi_lb = period_lb / 2 + delay_lb / (node_count * (node_count - 1))
    return (aoi=aoi_lb, period=period_lb, delay=delay_lb,
            independence_upper=independence_upper_bound,
            broadcasts_per_source=per_source_broadcast_lb)
end

function set_ieee_plot_font!()
    PyPlot.matplotlib.rcParams["font.family"] = "serif"
    PyPlot.matplotlib.rcParams["font.serif"] =
        ["Times New Roman", "Times", "DejaVu Serif"]
    PyPlot.matplotlib.rcParams["mathtext.fontset"] = "stix"
end

function plot_convergence(path, methods, incumbent_curves, current_curves,
                          theoretical_aoi_lb; systematic_aoi=nothing)
    set_ieee_plot_font!()
    figure(figsize=(7.2, 4.8))
    main_axis = gca()
    colors = Dict("Coverage-first" => "tab:blue", "Age-first" => "tab:orange",
                  "Dual-guided" => "tab:green")
    for method in methods
        main_axis.plot(1:length(current_curves[method]), current_curves[method],
                       color=colors[method], linestyle="--", alpha=0.28,
                       linewidth=0.9)
        main_axis.plot(1:length(incumbent_curves[method]), incumbent_curves[method],
                       color=colors[method], marker="o", markersize=3,
                       linewidth=1.6, label=method)
    end
    all_values = reduce(vcat, [vcat(incumbent_curves[method], current_curves[method])
                              for method in methods])
    main_axis.axhline(theoretical_aoi_lb, color="black", linestyle="--",
                      linewidth=1.3, label="topology-capacity LB")
    if systematic_aoi !== nothing
        main_axis.axhline(systematic_aoi, color="tab:purple", linestyle=":",
                          linewidth=1.5, label="sys (global interference)")
        main_axis.set_yscale("log")
        main_axis.set_ylim(0.88 * theoretical_aoi_lb, 1.15 * systematic_aoi)
    end
    xlabel("Multi-start / price-update iteration")
    ylabel("Best-so-far average AoI")
    PyPlot.grid(true, alpha=0.25)
    legend(fontsize=8, loc="center left", bbox_to_anchor=(0.01, 0.34))

    subplots_adjust(left=0.11, right=0.97, bottom=0.13, top=0.97)
    savefig(path, dpi=220)
    close()
end

function plot_decomposition(path, methods, best_results; systematic=nothing)
    set_ieee_plot_font!()
    figure(figsize=(7.2, 4.8))
    display_methods = systematic === nothing ? methods : [methods; "sys"]
    positions = 1:length(display_methods)
    periods = [best_results[method].period for method in methods]
    normalized_delays = [best_results[method].delay_sum /
                         (length(best_results[method].generation_times) *
                          (length(best_results[method].generation_times) - 1))
                         for method in methods]
    if systematic !== nothing
        push!(periods, systematic.period)
        push!(normalized_delays, systematic.mean_delay)
    end
    bar(positions .- 0.18, periods, width=0.36,
        label=L"update period $2\tau\;(=T)$")
    bar(positions .+ 0.18, normalized_delays, width=0.36,
        label=L"mean delay $\omega+\gamma$")
    xticks(positions, display_methods, rotation=18)
    systematic === nothing || yscale("log")
    ylabel("Slots")
    PyPlot.grid(true, axis="y", alpha=0.25)
    legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.01),
           ncol=2)
    subplots_adjust(left=0.11, right=0.97, bottom=0.20, top=0.86)
    savefig(path, dpi=220)
    close()
end

function plot_algorithm_workflow(path)
    set_ieee_plot_font!()
    figure(figsize=(12.5, 4.4))
    axis("off")
    function box(x, y, width, height, label; color="#EAF2F8")
        rectangle = PyPlot.matplotlib.patches.FancyBboxPatch(
            (x, y), width, height, boxstyle="round,pad=0.02",
            linewidth=1.5, edgecolor="#1F4E79", facecolor=color)
        gca().add_patch(rectangle)
        text(x + width / 2, y + height / 2, label, ha="center", va="center",
             fontsize=12, wrap=true)
    end
    function arrow(x1, y1, x2, y2; label="")
        annotate("", xy=(x2, y2), xytext=(x1, y1),
                 arrowprops=Dict("arrowstyle" => "->", "lw" => 1.5,
                                 "color" => "#34495E"))
        isempty(label) || text((x1+x2)/2, (y1+y2)/2 + 0.018, label,
                               ha="center", fontsize=10)
    end

    box(0.02, 0.59, 0.13, 0.16, "Network input\nC, source packets")
    box(0.18, 0.59, 0.17, 0.16, "Preprocessing\nGF and lower bounds")
    box(0.38, 0.59, 0.15, 0.16,
        "Period search\n\u0024T_{\\min},\\ldots,T_{\\max}\u0024")
    box(0.56, 0.54, 0.25, 0.26,
        "Exact fixed-T solver\nRMP  <->  source pricing\nclique cuts + branching",
        color="#FCF3CF")
    box(0.84, 0.59, 0.14, 0.16, "Output\nschedule, AoI, gap")

    box(0.35, 0.12, 0.27, 0.16,
        "Fast primal heuristics\nCoverage / Age / Dual\nincumbent UB + initial columns",
        color="#E8F8F5")

    arrow(0.15, 0.67, 0.18, 0.67)
    arrow(0.35, 0.67, 0.38, 0.67)
    arrow(0.53, 0.67, 0.56, 0.67, label="fixed T")
    arrow(0.81, 0.67, 0.84, 0.67, label="certified")
    arrow(0.58, 0.28, 0.64, 0.54)
    arrow(0.39, 0.28, 0.29, 0.59)
    text(0.615, 0.39, "UB + columns", ha="left", va="center", fontsize=10)
    text(0.335, 0.40, "warm start", ha="right", va="center", fontsize=10)

    xlim(0, 1)
    ylim(0.05, 0.86)
    title("Overall exact-and-anytime optimization architecture", fontsize=16)
    tight_layout()
    savefig(path, dpi=220, bbox_inches="tight")
    close()
end

function print_best_table(methods, best_results)
    println("\nBest result on the common topology")
    println(rpad("Method", 18), rpad("T", 8), rpad("Delay", 14),
            rpad("AoI", 14), rpad("Mean TX", 12), "Max TX")
    println(repeat("-", 76))
    for method in methods
        result = best_results[method]
        @printf("%-18s%-8d%-14d%-14.6f%-12.4f%d\n",
                method, result.period, result.delay_sum, result.average_aoi,
                result.mean_parallelism, result.maximum_parallelism)
    end
end

function comparison_main()
    node_count = env_int("AOI_N", 100)
    selected_source = env_int("AOI_SOURCE", 1)
    iterations = env_int("AOI_ITERATIONS", 20)
    seed = env_int("AOI_SEED", 2026)
    initial_radius = env_float("AOI_RADIUS", 0.18)
    methods = ["Coverage-first", "Age-first", "Dual-guided"]

    graph, positions, radius = connected_geometric_graph(
        node_count; seed=seed, initial_radius=initial_radius)
    conflict_graph = build_conflict_graph(graph)
    theoretical_bound = common_aoi_lower_bound(graph, conflict_graph)
    systematic = run_systematic_baseline(graph; exact=false)

    best_results = Dict{String,ComparisonResult}()
    incumbent_curves = Dict{String,Vector{Float64}}()
    current_curves = Dict{String,Vector{Float64}}()
    all_results = ComparisonResult[]
    total_runtime = @elapsed begin
        for method in methods
            best, results, incumbent_curve, current_curve = run_method(
                graph, conflict_graph, method, iterations, seed)
            best_results[method] = best
            incumbent_curves[method] = incumbent_curve
            current_curves[method] = current_curve
            append!(all_results, results)
            @printf("Completed %-16s: best AoI %.6f at iteration %d\n",
                    method, best.average_aoi, best.iteration)
        end
    end

    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    plot_directory = joinpath(project_root, "plots")
    mkpath(plot_directory)
    convergence_path = joinpath(plot_directory, "algorithm_convergence.png")
    decomposition_path = joinpath(plot_directory, "algorithm_decomposition.png")
    workflow_path = joinpath(plot_directory, "algorithm_workflow.png")
    result_path = joinpath(plot_directory, "algorithm_comparison.csv")
    best_method = argmin(method -> best_results[method].average_aoi, methods)
    best = best_results[best_method]
    source_path = joinpath(plot_directory,
                           "best_$(replace(lowercase(best_method), "-" => "_"))" *
                           "_source_$(selected_source).png")

    plot_convergence(convergence_path, methods, incumbent_curves,
                     current_curves, theoretical_bound.aoi;
                     systematic_aoi=systematic.average_aoi)
    plot_decomposition(decomposition_path, methods, best_results;
                       systematic=systematic)
    plot_algorithm_workflow(workflow_path)
    save_comparison_csv(result_path, all_results)
    for method in methods
        filename = replace(lowercase(method), "-" => "_") *
                   "_100_schedule.csv"
        save_schedule_csv(joinpath(plot_directory, filename),
                          best_results[method].schedule)
    end
    save_schedule_csv(joinpath(plot_directory, "systematic_100_schedule.csv"),
                      systematic.schedule)
    save_schedule_csv(joinpath(plot_directory, "best_100node_schedule.csv"),
                      best.schedule)
    plot_source_dissemination(graph, positions, selected_source,
                              best.first_parent, best.first_reception_slot,
                              source_path)

    println("\nCommon graph: N=$node_count, communication edges=$(ne(graph)), " *
            "conflict edges=$(ne(conflict_graph)), radius=$radius")
    println("Iterations per method: $iterations")
    println("Distance delay lower bound: $(theoretical_bound.delay)")
    println("Conflict-graph independence DSATUR upper bound: " *
            "$(theoretical_bound.independence_upper)")
    println("Broadcasts/source lower bound: $(theoretical_bound.broadcasts_per_source)")
    println("Common period lower bound: $(theoretical_bound.period)")
    @printf("Common theoretical AoI lower bound: %.6f\n", theoretical_bound.aoi)
    @printf("Global-interference SYS: T=%d, mean delay=%.6f, AoI=%.6f\n",
            systematic.period, systematic.mean_delay, systematic.average_aoi)
    print_best_table(methods, best_results)
    println("Best method: $best_method")
    @printf("Total comparison runtime: %.3f seconds\n", total_runtime)
    println("Convergence figure: $convergence_path")
    println("Decomposition figure: $decomposition_path")
    println("Algorithm workflow figure: $workflow_path")
    println("Comparison CSV: $result_path")
    println("Best source-path figure: $source_path")
end

function plot_existing_results_main()
    node_count = env_int("AOI_N", 100)
    seed = env_int("AOI_SEED", 2026)
    initial_radius = env_float("AOI_RADIUS", 0.18)
    methods = ["Coverage-first", "Age-first", "Dual-guided"]
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    plot_directory = joinpath(project_root, "plots")
    csv_path = joinpath(plot_directory, "algorithm_comparison.csv")
    lines = readlines(csv_path)[2:end]
    rows = [split(line, ',') for line in lines]
    incumbent_curves = Dict{String,Vector{Float64}}()
    current_curves = Dict{String,Vector{Float64}}()
    best_results = Dict{String,Any}()
    for method in methods
        method_rows = [row for row in rows if row[1] == method]
        sort!(method_rows, by=row -> parse(Int, row[2]))
        current = [parse(Float64, row[5]) for row in method_rows]
        incumbent = accumulate(min, current)
        incumbent_curves[method] = incumbent
        current_curves[method] = current
        best_index = argmin(current)
        best_row = method_rows[best_index]
        best_results[method] = (
            period=parse(Int, best_row[3]),
            delay_sum=parse(Int, best_row[4]),
            generation_times=zeros(Int, node_count))
    end
    graph, _, _ = connected_geometric_graph(
        node_count; seed=seed, initial_radius=initial_radius)
    conflict_graph = build_conflict_graph(graph)
    theoretical_bound = common_aoi_lower_bound(graph, conflict_graph)
    systematic = run_systematic_baseline(graph; exact=false)
    plot_convergence(joinpath(plot_directory, "algorithm_convergence.png"),
                     methods, incumbent_curves, current_curves,
                     theoretical_bound.aoi;
                     systematic_aoi=systematic.average_aoi)
    plot_decomposition(joinpath(plot_directory, "algorithm_decomposition.png"),
                       methods, best_results; systematic=systematic)
    println("Replotted existing 20-iteration CSV with LB=$(theoretical_bound.aoi)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if get(ENV, "AOI_WORKFLOW_ONLY", "0") == "1"
        project_root = normpath(joinpath(@__DIR__, "..", ".."))
        workflow_path = joinpath(project_root, "plots", "algorithm_workflow.png")
        plot_algorithm_workflow(workflow_path)
        println("Algorithm workflow figure: $workflow_path")
    elseif get(ENV, "AOI_PLOT_RESULTS_ONLY", "0") == "1"
        plot_existing_results_main()
    else
        comparison_main()
    end
end
