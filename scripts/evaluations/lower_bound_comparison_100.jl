include(joinpath(@__DIR__, "interference_algorithm_comparison.jl"))

using Statistics

function timed_value(function_to_run)
    elapsed = @elapsed value = function_to_run()
    return value, 1000.0 * elapsed
end

function distance_statistics(graph)
    node_count = nv(graph)
    distance_sum = 0
    diameter = 0
    for source in 1:node_count
        distances = gdistances(graph, source)
        for destination in 1:node_count
            source == destination && continue
            distance = distances[destination]
            distance_sum += distance
            diameter = max(diameter, distance)
        end
    end
    return distance_sum, diameter
end

function degree_sequence_broadcast_bound(graph)
    node_count = nv(graph)
    remaining = node_count - 1
    sorted_degrees = sort(degree(graph), rev=true)
    count = 0
    covered_upper = 0
    while covered_upper < remaining
        count += 1
        covered_upper += sorted_degrees[count]
    end
    return node_count * count, count
end

function source_launch_degree_sequence_bound(graph)
    node_count = nv(graph)
    degrees = degree(graph)
    total = 0
    per_source = Int[]
    for source in 1:node_count
        remaining = max(0, node_count - 1 - degrees[source])
        count = 1
        covered_upper = 0
        sorted_degrees = sort([degrees[node] for node in 1:node_count if node != source],
                              rev=true)
        index = 1
        while covered_upper < remaining
            covered_upper += sorted_degrees[index]
            count += 1
            index += 1
        end
        total += count
        push!(per_source, count)
    end
    return total, mean(per_source), maximum(per_source)
end

function best_heuristic_aoi(project_root)
    path = joinpath(project_root, "plots", "algorithm_comparison.csv")
    rows = readlines(path)[2:end]
    values = [parse(Float64, split(row, ',')[5]) for row in rows]
    return minimum(values)
end

function write_rows(path, rows)
    open(path, "w") do output
        println(output, "bound,period_lb,delay_lb,average_aoi_lb,compute_ms," *
                        "heuristic_to_lb,notes")
        for row in rows
            println(output, join(row, ','))
        end
    end
end

function main()
    project_root = normpath(joinpath(@__DIR__, "..", ".."))
    graph, _, radius = connected_geometric_graph(100; seed=2026, initial_radius=0.18)
    conflict = build_conflict_graph(graph)
    pair_count = nv(graph) * (nv(graph) - 1)
    best_aoi = best_heuristic_aoi(project_root)

    (distance_sum, diameter), distance_ms = timed_value(() -> distance_statistics(graph))
    distance_average = distance_sum / pair_count

    independence_upper, dsatur_ms =
        timed_value(() -> complement_dsatur_coloring_upper_bound(conflict))

    (common_bound, common_ms) = timed_value(() -> common_aoi_lower_bound(graph, conflict))

    (degree_total, degree_per_source), degree_ms =
        timed_value(() -> degree_sequence_broadcast_bound(graph))
    degree_period = ceil(Int, degree_total / independence_upper - 1e-9)
    degree_aoi = degree_period / 2 + distance_average

    (launch_total, launch_mean, launch_max), launch_ms =
        timed_value(() -> source_launch_degree_sequence_bound(graph))
    launch_period = ceil(Int, launch_total / independence_upper - 1e-9)
    launch_aoi = launch_period / 2 + distance_average

    hop_period = diameter
    hop_aoi = hop_period / 2 + distance_average

    rows = [
        ("distance-only", 0, distance_sum, distance_average, distance_ms,
         best_aoi / distance_average,
         "all-pairs BFS; no period-capacity term"),
        ("diameter-plus-distance", hop_period, distance_sum, hop_aoi,
         distance_ms, best_aoi / hop_aoi,
         "adds causal horizon diameter"),
        ("max-degree capacity", common_bound.period, distance_sum,
         common_bound.aoi, common_ms, best_aoi / common_bound.aoi,
         "current plotted topology-capacity LB"),
        ("degree-sequence capacity", degree_period, distance_sum, degree_aoi,
         distance_ms + dsatur_ms + degree_ms, best_aoi / degree_aoi,
         "uses sorted degrees instead of maximum degree"),
        ("launch-aware degree-sequence", launch_period, distance_sum, launch_aoi,
         distance_ms + dsatur_ms + launch_ms, best_aoi / launch_aoi,
         "accounts for mandatory source launch"),
    ]

    output_path = joinpath(project_root, "data", "exp_res",
                           "lower_bound_comparison_100.csv")
    write_rows(output_path, rows)

    @printf("100-node graph: radius=%.3f edges=%d conflict_edges=%d Δ=%d\n",
            radius, ne(graph), ne(conflict), maximum(degree(graph)))
    @printf("distance_sum=%d distance_avg=%.6f diameter=%d α_upper=%d\n",
            distance_sum, distance_average, diameter, independence_upper)
    @printf("best heuristic AoI=%.6f\n", best_aoi)
    println("Wrote $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
