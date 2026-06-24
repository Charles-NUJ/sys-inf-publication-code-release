using DrWatson
@quickactivate "age-optimal-multisource-flooding"
using Graphs
include(joinpath(srcdir(), "utilities", "state_utilities_src_inf_v3.jl"))
include(joinpath(srcdir(), "utilities", "aoi_utilities_src_inf_v3.jl"))

function pan_connectivity()
    graph = SimpleGraph(5)
    for edge in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
        add_edge!(graph, edge...)
    end
    return Matrix{Int8}(adjacency_matrix(graph))
end

function unit_payload(node_count, source)
    payload = zeros(Int8, node_count)
    payload[source] = 1
    return payload
end

function age_first_best_schedule()
    return [
        [(2, 2)],
        [(4, 2)],
        [(3, 3)],
        [(2, 3)],
        [(4, 4)],
        [(2, 4)],
        [(5, 5), (1, 1)],
        [(2, 1)],
        [(4, 5)],
        [(2, 5)],
        [(4, 1)],
    ]
end

function to_action_slot(decoders, payload_slot)
    if length(payload_slot) == 1
        tx, source = payload_slot[1]
        payload = unit_payload(size(decoders, 1), source)
        action = coded_payload_to_action(decoders, tx, payload)
        ismissing(action) && error("No state-matrix action for ($(tx),$(source))")
        return Dict("at" => "g", "a" => action)
    end
    actions = Tuple{Int64,Int64}[]
    for (tx, source) in payload_slot
        payload = unit_payload(size(decoders, 1), source)
        action = coded_payload_to_action(decoders, tx, payload)
        ismissing(action) && error("No state-matrix action for ($(tx),$(source))")
        push!(actions, action)
    end
    return Dict("at" => "m", "a" => Tuple(actions))
end

function main()
    connectivity = pan_connectivity()
    node_count = size(connectivity, 1)
    decoders = zeros(Int8, node_count, node_count, node_count)
    for node in 1:node_count
        decoders[node, 1, node] = 1
    end
    first_encodings = fill(Int8(127), node_count)
    decoding_levels = fill(Int8(127), node_count, node_count)
    content_purged = falses(node_count)
    action_schedule = Dict{String,Any}[]

    for (step, payload_slot) in enumerate(age_first_best_schedule())
        action_slot = to_action_slot(decoders, payload_slot)
        push!(action_schedule, action_slot)
        payload_slot_before = action_to_coded_payload_all_simu(decoders, action_slot)
        good, _ = apply_action_inf!(connectivity, decoders, action_slot)
        good || error("State update rejected at step $step: $action_slot")

        if payload_slot_before["pt"] == "g"
            tx = payload_slot_before["t"]
            payload = payload_slot_before["p"]
            first_encodings[tx] == 127 && payload[tx] == 1 &&
                (first_encodings[tx] = step)
        else
            for (tx, payload) in zip(payload_slot_before["t"], payload_slot_before["p"])
                first_encodings[tx] == 127 && payload[tx] == 1 &&
                    (first_encodings[tx] = step)
            end
        end

        decoded = get_is_decoded(decoders)
        for index in findall(decoded .== 1)
            decoding_levels[index] == 127 && (decoding_levels[index] = step + 1)
        end
        purge_disseminated_contents!(decoders, content_purged)
    end
    decoding_levels[I(node_count)] .= 1
    delay_sum = sum(decoding_levels[s, j] - first_encodings[s]
                    for s in 1:node_count, j in 1:node_count if s != j)
    aoi = length(age_first_best_schedule()) / 2 + delay_sum / (node_count * (node_count - 1))

    output_path = datadir("exp_res", "pan5_age_first_state_replay.txt")
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "state_action_schedule=$action_schedule")
        println(io, "first_encodings=$first_encodings")
        println(io, "first_decodings=$decoding_levels")
        println(io, "delay_sum=$delay_sum")
        println(io, "mean_delay=$(delay_sum / (node_count * (node_count - 1)))")
        println(io, "T=$(length(age_first_best_schedule()))")
        println(io, "aoi=$aoi")
    end

    println("state_action_schedule=$action_schedule")
    println("first_encodings=$first_encodings")
    println("first_decodings=$decoding_levels")
    println("delay_sum=$delay_sum")
    println("aoi=$aoi")
    println("Validation file: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
