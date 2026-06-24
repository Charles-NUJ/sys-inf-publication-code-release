using DrWatson
@quickactivate "age-optimal-multisource-flooding"
using Distributed
using Graphs
using LinearAlgebra
using Printf

include(joinpath(srcdir(), "solvers", "incremental_a_star_inf_v9.jl"))

function pan_connectivity()
    graph = SimpleGraph(5)
    for edge in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
        add_edge!(graph, edge...)
    end
    return Matrix{Int8}(adjacency_matrix(graph))
end

function format_payload_schedule(payload_schedule)
    slots = String[]
    for slot in payload_schedule
        if slot["pt"] == "g"
            source = findfirst(slot["p"] .== 1)
            push!(slots, "($(Int(slot["t"])),$source)")
        elseif slot["pt"] == "m"
            actions = String[]
            for (tx, payload) in zip(slot["t"], slot["p"])
                source = findfirst(payload .== 1)
                push!(actions, "($(Int(tx)),$source)")
            end
            push!(slots, "(" * join(actions, ",") * ")")
        else
            error("Unknown slot type $(slot["pt"])")
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

function write_result(path, solver_name, solution, runtime, max_inflight,
                      max_inference_simu, repair_starts)
    period, delay_sum, mean_delay, average_aoi, first_decodings =
        aoi_components(solution)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "solver=$solver_name")
        println(io, "max_inflight=$max_inflight")
        println(io, "max_inference_simu=$max_inference_simu")
        println(io, "repair_starts=$repair_starts")
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
    return period, delay_sum, mean_delay, average_aoi
end

function main()
    connectivity = pan_connectivity()
    solution_knowledge = Dict{String, Any}(
        "dissemination_delays" => zeros(5, 5),
        "t_star" => 12,
    )
    max_inflight = parse(Int, get(ENV, "AOI_PAN5_MAX_INFLIGHT", "2"))
    max_inference_simu = parse(Int, get(ENV, "AOI_PAN5_MAX_INFERENCE", "2"))
    repair_starts = parse(Int, get(ENV, "AOI_PAN5_REPAIR_STARTS", "500"))

    runtime_v8 = @elapsed solution_v8 = incremental_a_star_inf_v8_search(
        connectivity, 1, max_inflight, max_inference_simu,
        deepcopy(solution_knowledge);
        repair_starts = repair_starts)
    result_v8 = write_result(
        datadir("exp_res", "pan5_incremental_v8_validation.txt"),
        "incremental_a_star_inf_v8",
        solution_v8,
        runtime_v8,
        max_inflight,
        max_inference_simu,
        repair_starts)

    runtime_v9 = @elapsed solution_v9 = incremental_a_star_inf_v9_search(
        connectivity;
        max_coding_degree = 1,
        max_inflight = max_inflight,
        max_inference_simu = max_inference_simu,
        solution_knowledge = deepcopy(solution_knowledge),
        repair_starts = repair_starts)
    result_v9 = write_result(
        datadir("exp_res", "pan5_incremental_v9_validation.txt"),
        "incremental_a_star_inf_v9",
        solution_v9,
        runtime_v9,
        max_inflight,
        max_inference_simu,
        repair_starts)

    @printf("v8: T=%d delay_sum=%d mean_delay=%.4f AoI=%.4f runtime=%.3f\n",
            result_v8[1], result_v8[2], result_v8[3], result_v8[4], runtime_v8)
    println("v8 schedule=$(format_payload_schedule(solution_v8.payload_schedule))")
    @printf("v9: T=%d delay_sum=%d mean_delay=%.4f AoI=%.4f runtime=%.3f\n",
            result_v9[1], result_v9[2], result_v9[3], result_v9[4], runtime_v9)
    println("v9 schedule=$(format_payload_schedule(solution_v9.payload_schedule))")
    println("Validation files:")
    println(datadir("exp_res", "pan5_incremental_v8_validation.txt"))
    println(datadir("exp_res", "pan5_incremental_v9_validation.txt"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
