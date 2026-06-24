using DrWatson
@quickactivate "age-optimal-multisource-flooding"

using Graphs

function pan_connectivity()
    graph = SimpleGraph(5)
    for edge in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
        add_edge!(graph, edge...)
    end
    return Matrix{Int8}(adjacency_matrix(graph))
end

function unit_source(payload)
    nz = findall(payload .== 1)
    length(nz) == 1 || return missing
    return Int(nz[1])
end

function format_payload_schedule_global(payload_schedule)
    slots = String[]
    for (tx, payload) in payload_schedule
        source = unit_source(payload)
        ismissing(source) && error("Non-unit payload in global schedule: $payload")
        push!(slots, "($(Int(tx)),$source)")
    end
    return join(slots, ",")
end

function format_payload_schedule_inf(payload_schedule)
    slots = String[]
    for slot in payload_schedule
        if slot["pt"] == "g"
            source = unit_source(slot["p"])
            ismissing(source) && error("Non-unit payload in interference schedule: $(slot)")
            push!(slots, "($(Int(slot["t"])),$source)")
        elseif slot["pt"] == "m"
            actions = String[]
            for (tx, payload) in zip(slot["t"], slot["p"])
                source = unit_source(payload)
                ismissing(source) && error("Non-unit payload in interference schedule: $(slot)")
                push!(actions, "($(Int(tx)),$source)")
            end
            push!(slots, "(" * join(actions, ",") * ")")
        else
            error("Unknown slot type $(slot["pt"])")
        end
    end
    return join(slots, ",")
end

function aoi_from_events(first_encodings, decoding_levels, period)
    node_count = length(first_encodings)
    delay_sum = sum(decoding_levels[s, j] - first_encodings[s]
                    for s in 1:node_count, j in 1:node_count if s != j)
    return period, delay_sum, period / 2 + delay_sum / (node_count * (node_count - 1))
end

if abspath(PROGRAM_FILE) == @__FILE__
    connectivity = pan_connectivity()
    solution_knowledge = Dict{String,Any}(
        "dissemination_delays" => zeros(5, 5),
        "t_star" => 12,
    )

    include(joinpath(srcdir(), "solvers", "a_star_src.jl"))
    global_solution = a_star_search(connectivity, 1, 1, deepcopy(solution_knowledge))
    global_decoding = copy(global_solution.decoding_levels)
    global_decoding[I(5)] .= 1
    global_period, global_delay_sum, global_aoi =
        aoi_from_events(global_solution.first_encodings, global_decoding,
                        length(global_solution.schedule))

    include(joinpath(srcdir(), "solvers", "a_star_src_inf_v3.jl"))
    inf_solution = a_star_inf_search(connectivity, 1, 2, 2, deepcopy(solution_knowledge))
    inf_decoding = copy(inf_solution.decoding_levels)
    inf_decoding[I(5)] .= 1
    inf_period, inf_delay_sum, inf_aoi =
        aoi_from_events(inf_solution.first_encodings, inf_decoding,
                        length(inf_solution.schedule))

    output_path = datadir("exp_res", "pan5_astar_model_validation.txt")
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "global_solver=a_star_src.jl")
        println(io, "global_payload_schedule=$(format_payload_schedule_global(global_solution.payload_schedule))")
        println(io, "global_T=$global_period")
        println(io, "global_delay_sum=$global_delay_sum")
        println(io, "global_mean_delay=$(global_delay_sum / 20)")
        println(io, "global_aoi=$global_aoi")
        println(io, "global_first_encodings=$(global_solution.first_encodings)")
        println(io, "global_first_decodings=$global_decoding")
        println(io)
        println(io, "interference_solver=a_star_src_inf_v3.jl")
        println(io, "interference_payload_schedule=$(format_payload_schedule_inf(inf_solution.payload_schedule))")
        println(io, "interference_T=$inf_period")
        println(io, "interference_delay_sum=$inf_delay_sum")
        println(io, "interference_mean_delay=$(inf_delay_sum / 20)")
        println(io, "interference_aoi=$inf_aoi")
        println(io, "interference_first_encodings=$(inf_solution.first_encodings)")
        println(io, "interference_first_decodings=$inf_decoding")
    end

    println("Global A*: T=$global_period, delay_sum=$global_delay_sum, AoI=$global_aoi")
    println("Global schedule: $(format_payload_schedule_global(global_solution.payload_schedule))")
    println("Protocol-interference A*: T=$inf_period, delay_sum=$inf_delay_sum, AoI=$inf_aoi")
    println("Protocol schedule: $(format_payload_schedule_inf(inf_solution.payload_schedule))")
    println("Validation file: $output_path")
end
