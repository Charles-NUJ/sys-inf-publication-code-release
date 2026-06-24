using DrWatson
@quickactivate "age-optimal-multisource-flooding"

include(joinpath(@__DIR__, "interference_algorithm_comparison.jl"))

function env_int_local(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function validate_large_age_first()
    node_count = env_int_local("AOI_VALIDATE_N", 100)
    iterations = env_int_local("AOI_VALIDATE_ITERATIONS", 2)
    seed = env_int_local("AOI_VALIDATE_SEED", 2026)
    initial_radius = parse(Float64, get(ENV, "AOI_VALIDATE_RADIUS", "0.18"))

    graph, positions, radius = connected_geometric_graph(
        node_count; seed=seed, initial_radius=initial_radius)
    conflict = build_conflict_graph(graph)
    best, results, incumbent, current = run_method(
        graph, conflict, "Age-first", iterations, seed)
    theoretical_bound = common_aoi_lower_bound(graph, conflict)

    denominator = node_count * (node_count - 1)
    output_path = datadir("exp_res", "large_age_first_validation.csv")
    schedule_path = datadir("exp_res", "large_age_first_validation_schedule.csv")
    mkpath(dirname(output_path))

    open(output_path, "w") do io
        println(io, "nodes,iteration,period,delay_sum,mean_delay,average_aoi,recomputed_aoi,runtime,mean_parallelism,maximum_parallelism,lower_bound")
        for result in results
            mean_delay = result.delay_sum / denominator
            recomputed = result.period / 2 + mean_delay
            println(io, "$(node_count),$(result.iteration),$(result.period),$(result.delay_sum),$(mean_delay),$(result.average_aoi),$(recomputed),$(result.runtime),$(result.mean_parallelism),$(result.maximum_parallelism),$(theoretical_bound.aoi)")
        end
    end
    save_schedule_csv(schedule_path, best.schedule)

    mean_delay = best.delay_sum / denominator
    recomputed = best.period / 2 + mean_delay
    println("Age-first validation N=$node_count iterations=$iterations")
    println("edges=$(ne(graph)), conflict_edges=$(ne(conflict)), radius=$radius")
    println("best_iteration=$(best.iteration), T=$(best.period), delay_sum=$(best.delay_sum)")
    println("mean_delay=$mean_delay, average_aoi=$(best.average_aoi), recomputed=$recomputed")
    println("lower_bound=$(theoretical_bound.aoi)")
    println("validation_csv=$output_path")
    println("schedule_csv=$schedule_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_large_age_first()
end
