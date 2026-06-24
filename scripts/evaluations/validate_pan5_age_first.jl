using Graphs
using Random
using Printf

include(joinpath(@__DIR__, "interference_all_to_all_demo.jl"))

function pan_graph()
    graph = SimpleGraph(5)
    for (left, right) in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
        add_edge!(graph, left, right)
    end
    return graph
end

function format_schedule(schedule)
    slots = String[]
    for actions in schedule
        if length(actions) == 1
            action = actions[1]
            push!(slots, "($(action.transmitter),$(action.source))")
        else
            push!(slots, "(" * join(["($(action.transmitter),$(action.source))"
                                      for action in actions], ",") * ")")
        end
    end
    return join(slots, ",")
end

function delay_matrix(generation_times, arrival_times)
    node_count = length(generation_times)
    delays = zeros(Int, node_count, node_count)
    for source in 1:node_count, destination in 1:node_count
        if source != destination
            delays[source, destination] =
                arrival_times[source, destination] - generation_times[source]
        end
    end
    return delays
end

function run_age_first_once(graph, conflict, iteration)
    noise_scale = iteration == 1 ? 0.0 : 0.65
    return schedule_all_to_all(
        graph, conflict; policy=:age,
        rng=MersenneTwister(10_000 + iteration),
        noise_scale=noise_scale)
end

function main()
    graph = pan_graph()
    conflict = build_conflict_graph(graph)
    starts = parse(Int, get(ENV, "AOI_PAN5_AGE_STARTS", "50"))
    rows = NamedTuple[]
    best = nothing
    for iteration in 1:starts
        schedule, generation_times, arrival_times,
        first_parent, first_reception_slot =
            run_age_first_once(graph, conflict, iteration)
        validate_schedule(schedule, conflict)
        period, delay_sum, average_aoi =
            evaluate_schedule(schedule, generation_times, arrival_times)
        result = (iteration=iteration, schedule=schedule,
                  generation_times=generation_times,
                  arrival_times=arrival_times,
                  period=period, delay_sum=delay_sum,
                  average_aoi=average_aoi)
        push!(rows, result)
        if best === nothing || result.average_aoi < best.average_aoi
            best = result
        end
    end

    schedule = best.schedule
    generation_times = best.generation_times
    arrival_times = best.arrival_times
    period = best.period
    delay_sum = best.delay_sum
    average_aoi = best.average_aoi
    denominator = nv(graph) * (nv(graph) - 1)
    mean_delay = delay_sum / denominator

    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    output_path = joinpath(project_root, "data", "exp_res",
                           "pan5_age_first_validation.txt")
    iterations_path = joinpath(project_root, "data", "exp_res",
                               "pan5_age_first_iterations.csv")
    schedule_path = joinpath(project_root, "data", "exp_res",
                             "pan5_age_first_schedule.csv")
    mkpath(dirname(output_path))
    open(iterations_path, "w") do io
        println(io, "iteration,T,delay_sum,mean_delay,average_aoi,schedule")
        for row in rows
            row_mean_delay = row.delay_sum / denominator
            println(io, "$(row.iteration),$(row.period),$(row.delay_sum)," *
                        "$(row_mean_delay),$(row.average_aoi)," *
                        "\"$(format_schedule(row.schedule))\"")
        end
    end
    open(output_path, "w") do io
        println(io, "method=Age-first")
        println(io, "starts=$starts")
        println(io, "best_iteration=$(best.iteration)")
        println(io, "schedule=$(format_schedule(schedule))")
        println(io, "T=$period")
        println(io, "delay_sum=$delay_sum")
        println(io, "mean_delay=$mean_delay")
        println(io, "average_aoi=$average_aoi")
        println(io, "generation_times=$generation_times")
        println(io, "arrival_times=$arrival_times")
        println(io, "delay_matrix=$(delay_matrix(generation_times, arrival_times))")
        println(io, "previous_astar_log_aoi=7.45")
        println(io, "previous_astar_log_schedule=((1,1),(5,5)),(3,5),(2,1),(3,1),(2,2),(3,2),(2,5),(4,4),(2,4),(3,3),(2,3)")
    end
    save_schedule_csv(schedule_path, schedule)

    println("Age-first pan5 best of $starts starts")
    println("best_iteration=$(best.iteration)")
    println("schedule=$(format_schedule(schedule))")
    println("generation_times=$generation_times")
    println("arrival_times=$arrival_times")
    println("delay_matrix=$(delay_matrix(generation_times, arrival_times))")
    @printf("T=%d, delay_sum=%d, mean_delay=%.4f, AoI=%.4f\n",
            period, delay_sum, mean_delay, average_aoi)
    @printf("Difference to previous a_star_src_inf_v3 log: %.4f\n", average_aoi - 7.45)
    println("Validation file: $output_path")
    println("Iterations CSV: $iterations_path")
    println("Schedule CSV: $schedule_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
