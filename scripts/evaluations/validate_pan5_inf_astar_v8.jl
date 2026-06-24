using DrWatson
@quickactivate "age-optimal-multisource-flooding"
using Graphs
using LinearAlgebra
using Printf

include(joinpath(srcdir(), "solvers", "a_star_src_inf_v8.jl"))

function pan_connectivity()
    graph = SimpleGraph(5)
    for edge in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
        add_edge!(graph, edge...)
    end
    return Matrix{Int8}(adjacency_matrix(graph))
end

function payload_tuple(payload)
    nonzero = findall(payload .== 1)
    length(nonzero) == 1 || return missing
    return Int(nonzero[1])
end

function format_payload_schedule(payload_schedule)
    slots = String[]
    for slot in payload_schedule
        if slot["pt"] == "g"
            source = payload_tuple(slot["p"])
            ismissing(source) && error("Non-unit payload in slot: $slot")
            push!(slots, "($(Int(slot["t"])),$source)")
        elseif slot["pt"] == "m"
            actions = String[]
            for (tx, payload) in zip(slot["t"], slot["p"])
                source = payload_tuple(payload)
                ismissing(source) && error("Non-unit payload in slot: $slot")
                push!(actions, "($(Int(tx)),$source)")
            end
            push!(slots, "(" * join(actions, ",") * ")")
        else
            error("Unknown payload schedule slot type: $(slot["pt"])")
        end
    end
    return join(slots, ",")
end

function aoi_components(solution)
    node_count = length(solution.first_encodings)
    first_decodings = copy(solution.decoding_levels)
    first_decodings[I(node_count)] .= 1
    delay_sum = sum(first_decodings[source, destination] -
                    solution.first_encodings[source]
                    for source in 1:node_count,
                        destination in 1:node_count if source != destination)
    period = length(solution.payload_schedule)
    mean_delay = delay_sum / (node_count * (node_count - 1))
    average_aoi = period / 2 + mean_delay
    return period, delay_sum, mean_delay, average_aoi, first_decodings
end

function main()
    connectivity = pan_connectivity()
    node_count = size(connectivity, 1)
    solution_knowledge = Dict{String, Any}(
        "dissemination_delays" => zeros(5, 5),
        "t_star" => 12,
    )
    trace_path = datadir("exp_res", "astar_inf_v8_traces", "pan5_inf_astar_v8_traversal.csv")
    runtime = @elapsed solution = a_star_inf_v8_search(
        connectivity, 1, node_count, 2, deepcopy(solution_knowledge);
        save_traversal = true,
        trace_path = trace_path,
        trace_every = 1)
    period, delay_sum, mean_delay, average_aoi, first_decodings =
        aoi_components(solution)

    output_path = datadir("exp_res", "pan5_inf_astar_v8_validation.txt")
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "solver=a_star_src_inf_v8")
        println(io, "max_inflight=$node_count")
        println(io, "max_inference_simu=2")
        println(io, "trace_path=$trace_path")
        println(io, "payload_schedule=$(format_payload_schedule(solution.payload_schedule))")
        println(io, "raw_action_schedule=$(solution.schedule)")
        println(io, "first_encodings=$(solution.first_encodings)")
        println(io, "first_decodings=$first_decodings")
        println(io, "T=$period")
        println(io, "delay_sum=$delay_sum")
        println(io, "mean_delay=$mean_delay")
        println(io, "average_aoi=$average_aoi")
        println(io, "runtime_seconds=$runtime")
    end

    @printf("A* inf v8 pan5: T=%d delay_sum=%d mean_delay=%.4f AoI=%.4f runtime=%.3f\n",
            period, delay_sum, mean_delay, average_aoi, runtime)
    println("schedule=$(format_payload_schedule(solution.payload_schedule))")
    println("Validation file: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
