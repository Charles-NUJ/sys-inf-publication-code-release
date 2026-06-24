using DrWatson
@quickactivate "age-optimal-multisource-flooding"

include(joinpath(scriptsdir(), "evaluations", "top_src_inf_graph_analysis_inf.jl"))

function payload_tuple(payload)
    nonzero = findall(payload .== 1)
    length(nonzero) == 1 || return missing
    return Int(nonzero[1])
end

function payload_level_schedule(payload_schedule)
    slots = Vector{Vector{Tuple{Int,Int}}}()
    for slot in payload_schedule
        if slot["pt"] == "g"
            source = payload_tuple(slot["p"])
            ismissing(source) && error("Non-unit payload in validation: $(slot)")
            push!(slots, [(Int(slot["t"]), source)])
        elseif slot["pt"] == "m"
            actions = Tuple{Int,Int}[]
            for (tx, payload) in zip(slot["t"], slot["p"])
                source = payload_tuple(payload)
                ismissing(source) && error("Non-unit payload in validation: $(slot)")
                push!(actions, (Int(tx), source))
            end
            push!(slots, actions)
        else
            error("Unknown payload schedule slot type: $(slot)")
        end
    end
    return slots
end

function format_payload_schedule(slots)
    rendered = String[]
    for actions in slots
        if length(actions) == 1
            push!(rendered, "($(actions[1][1]),$(actions[1][2]))")
        else
            push!(rendered, "(" * join(["($(tx),$(source))" for (tx, source) in actions], ",") * ")")
        end
    end
    return join(rendered, ",")
end

function validate_pan5_inf_astar()
    config = Dict{String,Any}(
        "piggyback_enable" => false,
        "max_coding_degree" => 1,
        "graph_id" => 1,
        "max_inflight" => 2,
        "cutoff" => "",
        "eps" => 0.0,
        "graph_size" => 5,
        "heuristic" => false,
        "graph_type" => "pan",
        "inf" => true,
        "max_inference_simu" => 2,
        "regenerate_random" => false,
    )
    run_experiment(5, config)

    result_file = datadir("exp_raw_sys_inf", "pan", "5", savename(config, "jld2"))
    result = wload(result_file)
    slots = payload_level_schedule(result["payload_schedule"])
    denominator = 5 * 4
    delay_sum = sum(result["decoding_levels"][s, j] - result["first_encodings"][s]
                    for s in 1:5, j in 1:5 if s != j)
    half_period = result["t"] / 2
    mean_delay = delay_sum / denominator
    average_aoi = half_period + mean_delay

    output_path = datadir("exp_res", "pan5_inf_astar_validation.txt")
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        println(io, "payload_schedule=$(format_payload_schedule(slots))")
        println(io, "first_encodings=$(result["first_encodings"])")
        println(io, "decoding_levels=$(result["decoding_levels"])")
        println(io, "T=$(result["t"])")
        println(io, "delay_sum=$delay_sum")
        println(io, "mean_delay=$mean_delay")
        println(io, "average_aoi=$average_aoi")
        println(io, "aoi_get_new=$(result["aoi_get_new"])")
        println(io, "aoi_by_trapezoid=$(result["aoi_by_trapezoid"])")
        println(io, "result_file=$result_file")
    end

    println("Payload-level schedule: $(format_payload_schedule(slots))")
    println("T=$(result["t"]), delay_sum=$delay_sum, mean_delay=$mean_delay, AoI=$average_aoi")
    println("Result file: $result_file")
    println("Validation file: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_pan5_inf_astar()
end
