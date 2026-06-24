using Graphs
using Printf

"""One uncoded broadcast action: transmitter sends one source packet."""
struct BroadcastAction
    source::Int
    transmitter::Int
end

"""Create a connected N-node star with node 1 as the hub."""
function make_star_graph(node_count)
    graph = SimpleGraph(node_count)
    for leaf in 2:node_count
        add_edge!(graph, 1, leaf)
    end
    return graph
end

"""
Construct the certificate-optimal all-to-all schedule for a star.

The hub packet needs one hub broadcast. Every leaf packet needs exactly two
broadcasts: leaf -> hub, followed by hub -> all other leaves. Sources are
scheduled serially because generation time is the first source transmission;
serial shifting therefore adds no dissemination delay.
"""
function build_star_schedule(graph; hub=1)
    schedule = BroadcastAction[]
    for source in vertices(graph)
        if source == hub
            push!(schedule, BroadcastAction(source, hub))
        else
            push!(schedule, BroadcastAction(source, source))
            push!(schedule, BroadcastAction(source, hub))
        end
    end
    return schedule
end

"""Simulate broadcasts and return generation times, arrivals, and delay sum."""
function simulate_schedule(graph, schedule)
    node_count = nv(graph)
    informed = falses(node_count, node_count)
    generation_times = zeros(Int, node_count)
    arrival_times = zeros(Int, node_count, node_count)

    for source in 1:node_count
        informed[source, source] = true
    end

    for (time, action) in enumerate(schedule)
        source = action.source
        transmitter = action.transmitter
        informed[source, transmitter] ||
            error("Causality violation at t=$time: node $transmitter lacks source $source")

        if generation_times[source] == 0
            transmitter == source ||
                error("The first transmission of source $source is not made by its source")
            generation_times[source] = time
        end

        for receiver in neighbors(graph, transmitter)
            if !informed[source, receiver]
                informed[source, receiver] = true
                # Same convention as the original project: a slot-t packet is
                # decoded at timestamp t+1.
                arrival_times[source, receiver] = time + 1
            end
        end
    end

    all(informed) || error("The schedule does not complete all-to-all dissemination")
    delay_sum = sum(arrival_times[source, destination] - generation_times[source]
                    for source in 1:node_count, destination in 1:node_count
                    if source != destination)
    return generation_times, arrival_times, delay_sum
end

"""Universal lower bound: information needs at least graph distance hops."""
function distance_delay_lower_bound(graph)
    node_count = nv(graph)
    return sum(gdistances(graph, source)[destination]
               for source in 1:node_count, destination in 1:node_count
               if source != destination)
end

"""
Exact transmission lower bound for an uncoded star.

Each of N-1 leaf packets needs a leaf transmission and a hub relay. The hub
packet needs one hub transmission. No transmission can carry two sources.
"""
star_transmission_lower_bound(node_count) = 2 * (node_count - 1) + 1

function verify_star(graph; hub=1)
    node_count = nv(graph)
    degree(graph, hub) == node_count - 1 || return false
    return all(degree(graph, node) == 1 for node in vertices(graph) if node != hub)
end

function print_schedule_excerpt(schedule; count=6)
    println("First $count actions:")
    for time in 1:min(count, length(schedule))
        action = schedule[time]
        println("  t=$time: source=$(action.source), transmitter=$(action.transmitter)")
    end
    println("Last $count actions:")
    for time in max(1, length(schedule) - count + 1):length(schedule)
        action = schedule[time]
        println("  t=$time: source=$(action.source), transmitter=$(action.transmitter)")
    end
end

function main()
    node_count = 100
    hub = 1
    graph = make_star_graph(node_count)
    verify_star(graph; hub=hub) || error("Graph is not the expected star")

    elapsed = @elapsed begin
        schedule = build_star_schedule(graph; hub=hub)
        generation_times, arrival_times, delay_sum = simulate_schedule(graph, schedule)
        delay_lower_bound = distance_delay_lower_bound(graph)
        period_lower_bound = star_transmission_lower_bound(node_count)

        period = length(schedule)
        pair_count = node_count * (node_count - 1)
        average_aoi = period / 2 + delay_sum / pair_count
        aoi_lower_bound = period_lower_bound / 2 + delay_lower_bound / pair_count

        delay_optimal = delay_sum == delay_lower_bound
        period_optimal = period == period_lower_bound
        globally_optimal = delay_optimal && period_optimal && average_aoi == aoi_lower_bound

        println("Graph: connected star, N=$node_count, hub=$hub, edges=$(ne(graph))")
        println("All-to-all ordered source-destination pairs: $pair_count")
        println("Degree summary: hub=$(degree(graph, hub)), leaves=1")
        println("Top-left 10x10 block of the connectivity matrix:")
        display(Int.(Matrix(adjacency_matrix(graph)))[1:10, 1:10])
        println()
        println("Schedule period T: $period")
        println("Certified period lower bound: $period_lower_bound")
        println("Total first-packet delay: $delay_sum")
        println("All-pairs distance lower bound: $delay_lower_bound")
        @printf("Average AoI: %.10f\n", average_aoi)
        @printf("Certified AoI lower bound: %.10f\n", aoi_lower_bound)
        println("Period certificate attained: $period_optimal")
        println("Delay certificate attained: $delay_optimal")
        println("STRICT GLOBAL OPTIMUM CERTIFIED: $globally_optimal")
        println("Leaf-source delay example (source 2): " *
                "to hub=$(arrival_times[2, hub] - generation_times[2]), " *
                "to leaf 3=$(arrival_times[2, 3] - generation_times[2])")
        println("Hub-source delay example (source 1 to leaf 2): " *
                "$(arrival_times[1, 2] - generation_times[1])")
        print_schedule_excerpt(schedule)
    end
    @printf("Total construction + verification time: %.6f seconds\n", elapsed)
end

main()
