using Graphs
using Random
using Printf
using PyPlot

struct ScheduledBroadcast
    source::Int
    transmitter::Int
end

const DEFAULT_SCORE_PARAMETERS = (
    coverage_receiver=8.0,
    coverage_outstanding=0.25,
    coverage_launched=1.0,
    age_slope=0.03,
    launch_factor=0.20,
    age_receiver=2.0,
    dual_receiver=2.5,
)

function env_int(name, default)
    return parse(Int, get(ENV, name, string(default)))
end

function env_float(name, default)
    return parse(Float64, get(ENV, name, string(default)))
end

"""Generate a deterministic connected random geometric graph."""
function connected_geometric_graph(node_count; seed=2026, initial_radius=0.18,
                                   radius_step=0.01)
    Random.seed!(seed)
    positions = rand(2, node_count)
    radius = initial_radius

    while radius <= sqrt(2.0)
        graph = SimpleGraph(node_count)
        for left in 1:node_count-1, right in left+1:node_count
            distance = hypot(positions[1, left] - positions[1, right],
                             positions[2, left] - positions[2, right])
            distance <= radius && add_edge!(graph, left, right)
        end
        is_connected(graph) && return graph, positions, radius
        radius += radius_step
    end
    error("Failed to generate a connected geometric graph")
end

"""
Build the two-hop protocol-interference graph.

Two transmitters conflict if they are communication neighbors or if they have
a common communication neighbor. Thus feasible transmitters in one slot form
an independent set of this graph.
"""
function build_conflict_graph(communication_graph)
    node_count = nv(communication_graph)
    conflict_graph = SimpleGraph(node_count)
    neighbor_sets = [BitSet(neighbors(communication_graph, node))
                     for node in 1:node_count]

    for left in 1:node_count-1, right in left+1:node_count
        adjacent = has_edge(communication_graph, left, right)
        common_receiver = !isdisjoint(neighbor_sets[left], neighbor_sets[right])
        (adjacent || common_receiver) && add_edge!(conflict_graph, left, right)
    end
    return conflict_graph
end

function all_delivered(informed)
    return all(informed)
end

"""Return uninformed communication neighbors for a source/transmitter pair."""
function new_receivers(communication_graph, informed, source, transmitter)
    return [receiver for receiver in neighbors(communication_graph, transmitter)
            if !informed[source, receiver]]
end

"""
Choose at most one useful source packet per transmitter.

The score prioritizes packets with many outstanding destinations, broadcasts
that immediately cover many new receivers, and packets already in flight. The
result is a scalable primal heuristic and therefore an AoI upper bound.
"""
function transmitter_candidates(communication_graph, informed,
                                generation_times, time;
                                policy=:balanced,
                                rng=Random.default_rng(),
                                noise_scale=0.0,
                                node_prices=zeros(nv(communication_graph)),
                                score_parameters=DEFAULT_SCORE_PARAMETERS)
    node_count = nv(communication_graph)
    candidates = NamedTuple[]

    for transmitter in 1:node_count
        best = nothing
        best_score = -Inf
        for source in 1:node_count
            informed[source, transmitter] || continue
            receivers = new_receivers(communication_graph, informed,
                                      source, transmitter)
            isempty(receivers) && continue
            outstanding = node_count - count(informed[source, :])
            in_flight_age = generation_times[source] == 0 ? 0 :
                            time - generation_times[source]
            launched = generation_times[source] != 0
            if policy == :coverage
                score = score_parameters.coverage_receiver * length(receivers) +
                        score_parameters.coverage_outstanding * outstanding +
                        (launched ? score_parameters.coverage_launched : 0.0)
            elseif policy == :age
                age_pressure = launched ?
                    outstanding * (1.0 + score_parameters.age_slope * in_flight_age) :
                    score_parameters.launch_factor * outstanding
                score = age_pressure + score_parameters.age_receiver * length(receivers)
            elseif policy == :balanced
                launch_bonus = launched ? 2.0 : 0.0
                score = outstanding + 4.0 * length(receivers) +
                        0.05 * in_flight_age + launch_bonus
            elseif policy == :dual
                # Approximate the AoI stage cost (one unit per undelivered
                # destination and slot), then subtract a Lagrangian relay price.
                age_pressure = launched ?
                    outstanding * (1.0 + score_parameters.age_slope * in_flight_age) :
                    score_parameters.launch_factor * outstanding
                score = age_pressure + score_parameters.dual_receiver * length(receivers)
            else
                error("Unknown scheduling policy: $policy")
            end
            score -= node_prices[transmitter]
            score += noise_scale * randn(rng)
            if score > best_score
                best_score = score
                best = (source=source, transmitter=transmitter,
                        receivers=receivers, score=score)
            end
        end
        best === nothing || push!(candidates, best)
    end
    return candidates
end

"""Greedy weighted maximal independent set of candidate transmitters."""
function select_noninterfering(candidates, conflict_graph;
                               maximum_transmitters=typemax(Int))
    ordered = sort(candidates, by=candidate -> candidate.score, rev=true)
    selected = NamedTuple[]
    blocked = falses(nv(conflict_graph))

    for candidate in ordered
        transmitter = candidate.transmitter
        blocked[transmitter] && continue
        push!(selected, candidate)
        length(selected) >= maximum_transmitters && break
        blocked[transmitter] = true
        for conflict in neighbors(conflict_graph, transmitter)
            blocked[conflict] = true
        end
    end
    return selected
end

"""
Construct a feasible all-to-all schedule under two-hop interference.

All state updates in a slot are based on the state at the beginning of that
slot. The conflict graph guarantees that selected broadcasts have neither
adjacent transmitters nor a common receiver.
"""
function schedule_all_to_all(communication_graph, conflict_graph;
                             maximum_slots=1_000_000,
                             policy=:balanced,
                             rng=Random.default_rng(),
                             noise_scale=0.0,
                             node_prices=zeros(nv(communication_graph)),
                             score_parameters=DEFAULT_SCORE_PARAMETERS,
                             maximum_transmitters=typemax(Int))
    node_count = nv(communication_graph)
    informed = falses(node_count, node_count)
    generation_times = zeros(Int, node_count)
    arrival_times = zeros(Int, node_count, node_count)
    first_parent = zeros(Int, node_count, node_count)
    first_reception_slot = zeros(Int, node_count, node_count)
    schedule = Vector{Vector{ScheduledBroadcast}}()

    for source in 1:node_count
        informed[source, source] = true
    end

    time = 0
    while !all_delivered(informed)
        time += 1
        time > maximum_slots && error("Scheduler exceeded maximum_slots")
        candidates = transmitter_candidates(
            communication_graph, informed, generation_times, time;
            policy=policy, rng=rng, noise_scale=noise_scale,
            node_prices=node_prices, score_parameters=score_parameters)
        isempty(candidates) && error("No useful action exists; graph may be disconnected")
        selected = select_noninterfering(
            candidates, conflict_graph;
            maximum_transmitters=maximum_transmitters)
        isempty(selected) && error("Independent-set selector returned no action")

        slot_actions = ScheduledBroadcast[]
        # Record launches first, then apply all receptions simultaneously.
        for candidate in selected
            source = candidate.source
            transmitter = candidate.transmitter
            if generation_times[source] == 0
                transmitter == source ||
                    error("A packet was relayed before its source launch")
                generation_times[source] = time
            end
            push!(slot_actions, ScheduledBroadcast(source, transmitter))
        end

        for candidate in selected
            source = candidate.source
            transmitter = candidate.transmitter
            for receiver in candidate.receivers
                if !informed[source, receiver]
                    informed[source, receiver] = true
                    arrival_times[source, receiver] = time + 1
                    first_parent[source, receiver] = transmitter
                    first_reception_slot[source, receiver] = time
                end
            end
        end
        push!(schedule, slot_actions)
    end

    return schedule, generation_times, arrival_times,
           first_parent, first_reception_slot
end

function evaluate_schedule(schedule, generation_times, arrival_times)
    node_count = length(generation_times)
    period = length(schedule)
    delay_sum = sum(arrival_times[source, destination] - generation_times[source]
                    for source in 1:node_count, destination in 1:node_count
                    if source != destination)
    average_aoi = period / 2 + delay_sum / (node_count * (node_count - 1))
    return period, delay_sum, average_aoi
end

"""Universal delay lower bound based on communication-graph distances."""
function distance_delay_lower_bound(graph)
    node_count = nv(graph)
    return sum(gdistances(graph, source)[destination]
               for source in 1:node_count, destination in 1:node_count
               if source != destination)
end

function validate_schedule(schedule, conflict_graph)
    for (time, actions) in enumerate(schedule)
        transmitters = [action.transmitter for action in actions]
        length(unique(transmitters)) == length(transmitters) ||
            error("Node transmits two packets in slot $time")
        for left_index in 1:length(transmitters)-1
            for right_index in left_index+1:length(transmitters)
                has_edge(conflict_graph, transmitters[left_index],
                          transmitters[right_index]) &&
                    error("Interference violation in slot $time")
            end
        end
    end
    return true
end

function save_schedule_csv(path, schedule)
    open(path, "w") do output
        println(output, "slot,source,transmitter")
        for (time, actions) in enumerate(schedule), action in actions
            println(output, "$time,$(action.source),$(action.transmitter)")
        end
    end
end

function plot_source_dissemination(communication_graph, positions, source,
                                   first_parent, first_reception_slot, output_path)
    node_count = nv(communication_graph)
    reception_slots = first_reception_slot[source, :]
    maximum_reception = max(1, maximum(reception_slots))

    figure(figsize=(11, 9))
    for edge in edges(communication_graph)
        left, right = src(edge), dst(edge)
        plot([positions[1, left], positions[1, right]],
             [positions[2, left], positions[2, right]],
             color="0.85", linewidth=0.7, zorder=1)
    end

    labeled_events = Set{Tuple{Int,Int}}()
    for destination in 1:node_count
        destination == source && continue
        parent = first_parent[source, destination]
        parent == 0 && continue
        slot = first_reception_slot[source, destination]
        x1, y1 = positions[:, parent]
        x2, y2 = positions[:, destination]
        annotate("", xy=(x2, y2), xytext=(x1, y1),
                 arrowprops=Dict("arrowstyle" => "->", "color" => "tab:red",
                                 "lw" => 1.5, "alpha" => 0.8),
                 zorder=3)
        event = (parent, slot)
        if !(event in labeled_events)
            push!(labeled_events, event)
            text(x1 + 0.008, y1 - 0.012, "tx@$slot", fontsize=6,
                 color="darkred", zorder=4)
        end
    end

    colors = [node == source ? 0 : reception_slots[node]
              for node in 1:node_count]
    node_scatter = scatter(positions[1, :], positions[2, :], c=colors,
                           cmap="viridis", vmin=0, vmax=maximum_reception,
                           s=42, edgecolors="black", linewidths=0.4, zorder=5)
    scatter([positions[1, source]], [positions[2, source]], s=180,
            marker="*", color="gold", edgecolors="black", zorder=6,
            label="source $source")
    for node in 1:node_count
        text(positions[1, node] + 0.006, positions[2, node] + 0.006,
             string(node), fontsize=6, zorder=7)
    end

    colorbar(node_scatter, label="First-reception slot for selected source")
    legend(loc="best")
    xlabel("x")
    ylabel("y")
    axis("equal")
    PyPlot.grid(false)
    tight_layout()
    savefig(output_path, dpi=220)
    close()
end

function main()
    node_count = env_int("AOI_N", 100)
    selected_source = env_int("AOI_SOURCE", 1)
    seed = env_int("AOI_SEED", 2026)
    initial_radius = env_float("AOI_RADIUS", 0.18)
    1 <= selected_source <= node_count || error("AOI_SOURCE is outside 1:N")

    graph, positions, radius = connected_geometric_graph(
        node_count; seed=seed, initial_radius=initial_radius)
    conflict_graph = build_conflict_graph(graph)

    elapsed = @elapsed begin
        schedule, generation_times, arrival_times,
        first_parent, first_reception_slot = schedule_all_to_all(
            graph, conflict_graph)
        validate_schedule(schedule, conflict_graph)
        period, delay_sum, average_aoi = evaluate_schedule(
            schedule, generation_times, arrival_times)
        delay_lower_bound = distance_delay_lower_bound(graph)

        project_root = normpath(joinpath(@__DIR__, "..", ".."))
        plot_directory = joinpath(project_root, "plots")
        mkpath(plot_directory)
        figure_path = joinpath(plot_directory,
                               "interference_source_$(selected_source).png")
        csv_path = joinpath(plot_directory, "interference_schedule.csv")
        plot_source_dissemination(graph, positions, selected_source,
                                  first_parent, first_reception_slot, figure_path)
        save_schedule_csv(csv_path, schedule)

        simultaneous = [length(actions) for actions in schedule]
        @printf("Connected geometric graph: N=%d, edges=%d, radius=%.3f\n",
                node_count, ne(graph), radius)
        @printf("Two-hop conflict graph: edges=%d, density=%.5f\n",
                ne(conflict_graph), 2ne(conflict_graph) /
                (node_count * (node_count - 1)))
        println("Period T: $period")
        println("Total first-packet delay: $delay_sum")
        println("Distance-only delay lower bound: $delay_lower_bound")
        @printf("Delay ratio to distance lower bound: %.6f\n",
                delay_sum / delay_lower_bound)
        @printf("Average AoI: %.10f\n", average_aoi)
        @printf("Mean simultaneous transmitters: %.4f\n", sum(simultaneous) / period)
        println("Maximum simultaneous transmitters: $(maximum(simultaneous))")
        println("Selected source: $selected_source")
        println("Source generation slot: $(generation_times[selected_source])")
        println("Source completion timestamp: $(maximum(arrival_times[selected_source, :]))")
        println("Figure: $figure_path")
        println("Schedule CSV: $csv_path")
        println("Optimality note: this scalable schedule is a feasible upper bound; " *
                "strict optimality requires the paper's branch-price-and-cut gap to close.")
    end
    @printf("Scheduling + evaluation + plotting time: %.6f seconds\n", elapsed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
