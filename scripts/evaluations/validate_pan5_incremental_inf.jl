using DrWatson
@quickactivate "age-optimal-multisource-flooding"
using Graphs
using Distributed
include(joinpath(srcdir(), "solvers", "incremental_a_star_inf_v3.jl"))

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

function main()
    connectivity = pan_connectivity()
    solution_knowledge = Dict{String,Any}(
        "dissemination_delays" => zeros(5, 5),
        "t_star" => 12,
    )
    runtime = @elapsed solution = incremental_a_star_inf_search(
        connectivity, 1, 2, 2, deepcopy(solution_knowledge))
    first_decodings = copy(solution.decoding_levels)
    first_decodings[I(5)] .= 1
    delay_sum = sum(first_decodings[s, j] - solution.first_encodings[s]
                    for s in 1:5, j in 1:5 if s != j)
    period = length(solution.payload_schedule)
    mean_delay = delay_sum / 20
    aoi = period / 2 + mean_delay

    output_path = datadir("exp_res", "pan5_incremental_inf_validation.txt")
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "solver=incremental_a_star_inf_v3")
        println(io, "max_inflight=2")
        println(io, "payload_schedule=$(format_payload_schedule(solution.payload_schedule))")
        println(io, "raw_action_schedule=$(solution.schedule)")
        println(io, "first_encodings=$(solution.first_encodings)")
        println(io, "first_decodings=$first_decodings")
        println(io, "T=$period")
        println(io, "delay_sum=$delay_sum")
        println(io, "mean_delay=$mean_delay")
        println(io, "average_aoi=$aoi")
        println(io, "runtime_seconds=$runtime")
        println(io, "age_first_feasible_reference=7.3")
    end

    println("Incremental A* pan5")
    println("payload_schedule=$(format_payload_schedule(solution.payload_schedule))")
    println("T=$period delay_sum=$delay_sum mean_delay=$mean_delay AoI=$aoi runtime=$runtime")
    println("Validation file: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
