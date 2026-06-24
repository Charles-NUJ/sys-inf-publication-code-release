include(joinpath(@__DIR__, "interference_algorithm_comparison.jl"))

function scalable_radius(node_count)
    return 0.18 * sqrt(100 / node_count)
end

function run_age_first_scale(node_count, iterations, seed)
    initial_radius = scalable_radius(node_count)
    graph, positions, radius = connected_geometric_graph(
        node_count; seed=seed, initial_radius=initial_radius,
        radius_step=0.0025)
    conflict = build_conflict_graph(graph)
    best, results, incumbent, current = run_method(
        graph, conflict, "Age-first", iterations, seed)
    return (node_count=node_count, graph=graph, positions=positions,
            conflict=conflict, radius=radius, best=best, results=results,
            incumbent=incumbent, current=current)
end

function save_large_scale_data(output_directory, experiments)
    iteration_path = joinpath(output_directory,
                              "age_first_large_iterations.csv")
    summary_path = joinpath(output_directory, "age_first_large_summary.csv")
    open(iteration_path, "w") do output
        println(output, "nodes,iteration,period,mean_delay,average_aoi," *
                        "incumbent_aoi,runtime_seconds,mean_parallelism," *
                        "maximum_parallelism")
        for experiment in experiments
            denominator = experiment.node_count * (experiment.node_count - 1)
            for (index, result) in enumerate(experiment.results)
                mean_delay = result.delay_sum / denominator
                println(output, "$(experiment.node_count),$(result.iteration)," *
                    "$(result.period),$mean_delay,$(result.average_aoi)," *
                    "$(experiment.incumbent[index]),$(result.runtime)," *
                    "$(result.mean_parallelism),$(result.maximum_parallelism)")
            end
        end
    end
    open(summary_path, "w") do output
        println(output, "nodes,communication_edges,conflict_edges,radius," *
                        "best_iteration,period,mean_delay,average_aoi," *
                        "total_runtime_seconds,mean_parallelism," *
                        "maximum_parallelism")
        for experiment in experiments
            best = experiment.best
            denominator = experiment.node_count * (experiment.node_count - 1)
            mean_delay = best.delay_sum / denominator
            total_runtime = sum(result.runtime for result in experiment.results)
            println(output, "$(experiment.node_count),$(ne(experiment.graph))," *
                "$(ne(experiment.conflict)),$(experiment.radius)," *
                "$(best.iteration),$(best.period),$mean_delay," *
                "$(best.average_aoi),$total_runtime,$(best.mean_parallelism)," *
                "$(best.maximum_parallelism)")
        end
    end
    return iteration_path, summary_path
end

function plot_large_convergence(path, experiments)
    set_ieee_plot_font!()
    figure(figsize=(7.2, 4.6))
    colors = ["tab:purple", "tab:blue", "tab:orange", "tab:green"]
    for (index, experiment) in enumerate(experiments)
        first_aoi = experiment.incumbent[1]
        improvement = 100 .* (first_aoi .- experiment.incumbent) ./ first_aoi
        plot(1:length(improvement), improvement, marker="o", markersize=3.5,
             color=colors[index], linewidth=1.7,
             label="N=$(experiment.node_count)")
    end
    xlabel("Multi-start iteration")
    ylabel("Incumbent improvement from first start (%)")
    title("Age-first multi-start convergence")
    PyPlot.grid(true, alpha=0.25)
    legend(fontsize=8, loc="best")
    tight_layout()
    savefig(path, dpi=240)
    close()
end

function replot_existing_large_main()
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    data_directory = joinpath(project_root, "data", "exp_res")
    plot_directory = joinpath(project_root, "plots")
    iteration_rows = [split(line, ',') for line in
                      readlines(joinpath(data_directory,
                                         "age_first_large_iterations.csv"))[2:end]]
    summary_rows = [split(line, ',') for line in
                    readlines(joinpath(data_directory,
                                       "age_first_large_summary.csv"))[2:end]]
    experiments = Any[]
    comparison_rows = [split(line, ',') for line in
                       readlines(joinpath(project_root, "plots",
                                          "algorithm_comparison.csv"))[2:end]]
    age_rows = [row for row in comparison_rows if row[1] == "Age-first"]
    sort!(age_rows, by=row -> parse(Int, row[2]))
    current100 = [parse(Float64, row[5]) for row in age_rows]
    incumbent100 = accumulate(min, current100)
    best100_row = age_rows[argmin(current100)]
    graph100, _, radius100 = connected_geometric_graph(
        100; seed=2026, initial_radius=0.18)
    conflict100 = build_conflict_graph(graph100)
    push!(experiments, (
        node_count=100, current=current100, incumbent=incumbent100,
        best=(average_aoi=parse(Float64, best100_row[5]),
              period=parse(Int, best100_row[3]),
              delay_sum=parse(Float64, best100_row[4])),
        results=[(runtime=parse(Float64, row[6]),) for row in age_rows],
        communication_edges=ne(graph100), conflict_edges=ne(conflict100),
        radius=radius100, best_iteration=parse(Int, best100_row[2]),
        starts=length(age_rows)))
    for summary in summary_rows
        node_count = parse(Int, summary[1])
        rows = [row for row in iteration_rows if parse(Int, row[1]) == node_count]
        sort!(rows, by=row -> parse(Int, row[2]))
        current = [parse(Float64, row[5]) for row in rows]
        incumbent = [parse(Float64, row[6]) for row in rows]
        mean_delay = parse(Float64, summary[7])
        best = (average_aoi=parse(Float64, summary[8]),
                period=parse(Int, summary[6]),
                delay_sum=mean_delay * node_count * (node_count - 1))
        results = [(runtime=parse(Float64, summary[9]),)]
        push!(experiments, (node_count=node_count, current=current,
                            incumbent=incumbent, best=best, results=results,
                            communication_edges=parse(Int, summary[2]),
                            conflict_edges=parse(Int, summary[3]),
                            radius=parse(Float64, summary[4]),
                            best_iteration=parse(Int, summary[5]),
                            starts=length(rows)))
    end
    sort!(experiments, by=experiment -> experiment.node_count)
    combined_summary = joinpath(data_directory,
                                "age_first_100_1000_summary.csv")
    open(combined_summary, "w") do output
        println(output, "nodes,communication_edges,conflict_edges,radius,starts," *
                        "best_iteration,period,mean_delay,average_aoi," *
                        "total_runtime_seconds")
        for experiment in experiments
            denominator = experiment.node_count * (experiment.node_count - 1)
            mean_delay = experiment.best.delay_sum / denominator
            total_runtime = sum(result.runtime for result in experiment.results)
            println(output, "$(experiment.node_count)," *
                "$(experiment.communication_edges),$(experiment.conflict_edges)," *
                "$(experiment.radius),$(experiment.starts)," *
                "$(experiment.best_iteration),$(experiment.best.period)," *
                "$mean_delay,$(experiment.best.average_aoi),$total_runtime")
        end
    end
    plot_large_convergence(joinpath(plot_directory,
                                    "age_first_large_convergence.png"),
                           experiments)
    plot_large_results(joinpath(plot_directory, "age_first_large_results.png"),
                       experiments)
end

function plot_large_results(path, experiments)
    set_ieee_plot_font!()
    nodes = [experiment.node_count for experiment in experiments]
    aoi = [experiment.best.average_aoi for experiment in experiments]
    periods = [experiment.best.period for experiment in experiments]
    delays = [experiment.best.delay_sum /
              (experiment.node_count * (experiment.node_count - 1))
              for experiment in experiments]
    runtimes = [sum(result.runtime for result in experiment.results)
                for experiment in experiments]

    figure(figsize=(7.2, 4.8))
    left = gca()
    left.plot(nodes, aoi, marker="o", linewidth=1.8, label="average AoI")
    left.plot(nodes, [period / 2 for period in periods], marker="s",
              linewidth=1.5, label=L"half-period $T/2$")
    left.plot(nodes, delays, marker="^", linewidth=1.5,
              label="mean delay")
    left.set_xlabel("Number of nodes N")
    left.set_ylabel("Slots")
    left.grid(true, alpha=0.25)
    right = left.twinx()
    right.plot(nodes, runtimes, marker="D", color="tab:red", linewidth=1.8,
               linestyle="--", label="multi-start runtime")
    right.set_yscale("log")
    right.set_ylabel("Scheduling runtime (s, log scale)", color="tab:red")
    right.tick_params(axis="y", labelcolor="tab:red")
    left_lines, left_labels = left.get_legend_handles_labels()
    right_lines, right_labels = right.get_legend_handles_labels()
    left.legend(vcat(left_lines, right_lines), vcat(left_labels, right_labels),
                fontsize=8, loc="upper left")
    tight_layout()
    savefig(path, dpi=240)
    close()
end

function large_scale_main()
    iterations = env_int("AOI_ITERATIONS", 3)
    seed = env_int("AOI_SEED", 2026)
    node_counts = parse.(Int, split(get(ENV, "AOI_NODE_COUNTS", "400,700,1000"), ','))
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    output_directory = joinpath(project_root, "data", "exp_res")
    plot_directory = joinpath(project_root, "plots")
    mkpath(output_directory)
    mkpath(plot_directory)

    experiments = Any[]
    for node_count in node_counts
        @printf("Running Age-first: N=%d, iterations=%d\n", node_count, iterations)
        experiment = run_age_first_scale(node_count, iterations,
                                         seed + node_count)
        push!(experiments, experiment)
        @printf("N=%d: edges=%d, conflict=%d, radius=%.4f, best AoI=%.6f, T=%d\n",
                node_count, ne(experiment.graph), ne(experiment.conflict),
                experiment.radius, experiment.best.average_aoi,
                experiment.best.period)
    end

    iteration_path, summary_path = save_large_scale_data(
        output_directory, experiments)
    convergence_path = joinpath(plot_directory,
                                "age_first_large_convergence.png")
    results_path = joinpath(plot_directory, "age_first_large_results.png")
    plot_large_convergence(convergence_path, experiments)
    plot_large_results(results_path, experiments)
    println("Iterations: $iteration_path")
    println("Summary: $summary_path")
    println("Convergence: $convergence_path")
    println("Results: $results_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if get(ENV, "AOI_PLOT_RESULTS_ONLY", "0") == "1"
        replot_existing_large_main()
    else
        large_scale_main()
    end
end
