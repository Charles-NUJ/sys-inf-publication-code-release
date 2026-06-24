using Dates
using DataStructures

include(joinpath(@__DIR__, "a_star_src_inf_v3.jl"))

const ASTAR_INF_V8_TRACE_DIR = "astar_inf_v8_traces"

function _astar_inf_v8_make_simu_tx_list(connectivity, max_inference_simu)
	N = size(connectivity, 1)
	simu_tx_list = Vector{Vector{Int}}()
	if max_inference_simu <= 1
		return simu_tx_list
	end

	for node_mask in 1:(2^UInt32(N)-1)
		node_list = get_all_index(node_mask)
		if length(node_list) <= 1 || length(node_list) > max_inference_simu
			continue
		end

		feasible = true
		for first_index in 1:(length(node_list)-1)
			first_neighbors = findall(connectivity[node_list[first_index], :] .> 0)
			for second_index in (first_index+1):length(node_list)
				second_neighbors = findall(connectivity[node_list[second_index], :] .> 0)
				has_common_receiver = !isempty(intersect(first_neighbors, second_neighbors))
				adjacent_tx = node_list[first_index] in second_neighbors || node_list[second_index] in first_neighbors
				if has_common_receiver || adjacent_tx
					feasible = false
					break
				end
			end
			if !feasible
				break
			end
		end

		if feasible
			push!(simu_tx_list, node_list)
		end
	end
	return simu_tx_list
end

function _astar_inf_v8_start_state(connectivity)
	N = size(connectivity, 1)
	start_decoders = zeros(Int8, N, N, N)
	for node_id in 1:N
		start_decoders[node_id, 1, node_id] = 1
	end
	decoding_levels = fill(Int8(127), N, N)
	for node_id in 1:N
		decoding_levels[node_id, node_id] = 1
	end
	return AStarInfState(
		Dict{String, Any}[],
		Dict{String, Any}[],
		start_decoders,
		decoding_levels,
		fill(Int8(127), N),
		falses(N),
		0.0f0,
		0.0f0,
		0.0f0,
		0.0f0,
		false,
		Int8[],
	)
end

function _astar_inf_v8_inactive_nodes(state, true_leaves)
	inactive_nodes = Int[]
	for scheduled_action in state.schedule
		if scheduled_action["at"] == "g"
			action = scheduled_action["a"]
			union!(inactive_nodes, intersect(true_leaves, [action[1]]))
		elseif scheduled_action["at"] == "m"
			for action in scheduled_action["a"]
				union!(inactive_nodes, intersect(true_leaves, [action[1]]))
			end
		end
	end
	return inactive_nodes
end

function _astar_inf_v8_delay_sum(state)
	N = length(state.first_encodings)
	delay_sum = 0
	for source_id in 1:N
		generation_slot = state.first_encodings[source_id]
		if generation_slot == 127
			continue
		end
		for receiver_id in 1:N
			if source_id == receiver_id
				continue
			end
			receive_slot = state.decoding_levels[source_id, receiver_id]
			if receive_slot != 127
				delay_sum += Int(receive_slot) - Int(generation_slot)
			end
		end
	end
	return delay_sum
end

function _astar_inf_v8_scaled_lower_bound(state)
	N = length(state.first_encodings)
	normalizer = N * (N - 1)
	step = length(state.schedule)
	fixed_delay_sum = _astar_inf_v8_delay_sum(state)
	return step * normalizer / 2 + fixed_delay_sum
end

function _astar_inf_v8_scaled_exact_cost(state)
	if !state.end_state
		return Inf
	end
	return _astar_inf_v8_scaled_lower_bound(state)
end

function _astar_inf_v8_aoi(state)
	N = length(state.first_encodings)
	normalizer = N * (N - 1)
	return _astar_inf_v8_scaled_exact_cost(state) / normalizer
end

function _astar_inf_v8_decoded_pairs(state)
	N = length(state.first_encodings)
	decoded_pairs = 0
	for source_id in 1:N
		for receiver_id in 1:N
			if source_id != receiver_id && state.decoding_levels[source_id, receiver_id] != 127
				decoded_pairs += 1
			end
		end
	end
	return decoded_pairs
end

function _astar_inf_v8_action_to_string(action)
	if action === nothing
		return ""
	end
	return replace(string(action), "\"" => "\"\"")
end

function _astar_inf_v8_trace_path(N)
	trace_dir = datadir("exp_res", ASTAR_INF_V8_TRACE_DIR)
	mkpath(trace_dir)
	timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
	return joinpath(trace_dir, "astar_inf_v8_N$(N)_$(timestamp).csv")
end

function _astar_inf_v8_open_trace(save_traversal, trace_path, N)
	if !save_traversal
		return nothing, nothing
	end
	path = isnothing(trace_path) ? _astar_inf_v8_trace_path(N) : trace_path
	mkpath(dirname(path))
	io = open(path, "w")
	println(io, "event,iteration,state_id,parent_id,depth,lower_bound,exact_aoi,open_size,closed_size,decoded_pairs,delay_sum,action")
	return io, path
end

function _astar_inf_v8_write_trace!(io, event, iteration, state_id, parent_id, state, lower_bound, exact_aoi, open_size, closed_size, action)
	if isnothing(io)
		return nothing
	end
	action_text = _astar_inf_v8_action_to_string(action)
	println(io, join((
		event,
		iteration,
		state_id,
		parent_id,
		length(state.schedule),
		lower_bound,
		exact_aoi,
		open_size,
		closed_size,
		_astar_inf_v8_decoded_pairs(state),
		_astar_inf_v8_delay_sum(state),
		"\"$(action_text)\"",
	), ","))
	return nothing
end

function _astar_inf_v8_maybe_push_snapshot!(snapshots, state_id, parent_id, state, lower_bound, action)
	if isnothing(snapshots)
		return nothing
	end
	push!(snapshots, Dict(
		"state_id" => state_id,
		"parent_id" => parent_id,
		"depth" => length(state.schedule),
		"lower_bound" => lower_bound,
		"action" => deepcopy(action),
		"state" => deepcopy(state),
	))
	return nothing
end

function _astar_inf_v8_save_snapshots(snapshots, trace_path)
	if isnothing(snapshots) || isnothing(trace_path)
		return nothing
	end
	snapshot_path = replace(trace_path, r"\.csv$" => ".jld2")
	wsave(snapshot_path, Dict("snapshots" => snapshots))
	return snapshot_path
end

function _astar_inf_v8_schedule_payload_summary!(state)
	state.payload_schedule = Dict{String, Any}[]
	for scheduled_action in state.schedule
		payload = action_to_coded_payload_all_simu(state.decoders, scheduled_action)
		push!(state.payload_schedule, payload)
	end
	return state
end

function a_star_inf_v8_search(
	connectivity,
	max_coding_degree,
	max_inflight,
	max_inference_simu,
	solution_knowledge,
	star_epsilon = 1.0;
	save_traversal = false,
	trace_path = nothing,
	trace_every = 1,
	save_state_snapshots = false,
	max_expansions = typemax(Int),
)
	@info "no-e inf v8 exact traversal"
	_ = star_epsilon
	N = size(connectivity, 1)
	effective_max_inflight = isnothing(max_inflight) ? N : max_inflight
	_, _, true_leaves, _, _ = get_node_types(connectivity)
	action_set = vec(collect(Iterators.product(1:N, 1:2^N)))
	simu_tx_list = _astar_inf_v8_make_simu_tx_list(connectivity, max_inference_simu)

	start = _astar_inf_v8_start_state(connectivity)
	start_lower_bound = _astar_inf_v8_scaled_lower_bound(start)
	start.f_score = Float32(start_lower_bound)
	start.g_score = 0.0f0

	open_set_pq = PriorityQueue{Int, Float64}()
	state_by_id = Dict{Int, AStarInfState}()
	state_key_by_id = Dict{Int, UInt}()
	parent_by_id = Dict{Int, Int}()
	action_by_id = Dict{Int, Any}()
	best_lower_by_key = Dict{UInt, Float64}()
	closed_lower_by_key = Dict{UInt, Float64}()

	next_state_id = 1
	state_by_id[next_state_id] = start
	state_key_by_id[next_state_id] = hash(start)
	parent_by_id[next_state_id] = 0
	action_by_id[next_state_id] = nothing
	best_lower_by_key[hash(start)] = start_lower_bound
	open_set_pq[next_state_id] = start_lower_bound

	trace_io, actual_trace_path = _astar_inf_v8_open_trace(save_traversal, trace_path, N)
	snapshots = save_state_snapshots ? Vector{Dict{String, Any}}() : nothing
	_astar_inf_v8_write_trace!(trace_io, "start", 0, next_state_id, 0, start, start_lower_bound, "", length(open_set_pq), 0, nothing)
	_astar_inf_v8_maybe_push_snapshot!(snapshots, next_state_id, 0, start, start_lower_bound, nothing)

	best_solution = missing
	best_scaled_cost = Inf
	iteration = 0
	try
		while !isempty(open_set_pq)
			current_id = dequeue!(open_set_pq)
			current_state = state_by_id[current_id]
			current_key = state_key_by_id[current_id]
			current_lower_bound = _astar_inf_v8_scaled_lower_bound(current_state)

			if current_lower_bound >= best_scaled_cost
				break
			end
			if current_lower_bound > get(best_lower_by_key, current_key, Inf)
				continue
			end
			if current_lower_bound >= get(closed_lower_by_key, current_key, Inf)
				continue
			end

			iteration += 1
			closed_lower_by_key[current_key] = current_lower_bound
			if iteration % trace_every == 0
				_astar_inf_v8_write_trace!(trace_io, "expand", iteration, current_id, parent_by_id[current_id], current_state, current_lower_bound, "", length(open_set_pq), length(closed_lower_by_key), action_by_id[current_id])
			end

			if iteration > max_expansions
				@warn "A* v8 reached max_expansions before proving optimality" max_expansions best_scaled_cost current_lower_bound
				break
			end

			if current_state.end_state
				current_exact_cost = _astar_inf_v8_scaled_exact_cost(current_state)
				if current_exact_cost < best_scaled_cost
					best_scaled_cost = current_exact_cost
					best_solution = current_state
					_astar_inf_v8_write_trace!(trace_io, "incumbent", iteration, current_id, parent_by_id[current_id], current_state, current_lower_bound, _astar_inf_v8_aoi(current_state), length(open_set_pq), length(closed_lower_by_key), action_by_id[current_id])
				end
				continue
			end

			inactive_nodes = _astar_inf_v8_inactive_nodes(current_state, true_leaves)
			branching_actions = valid_actions_all_simu(
				action_set,
				simu_tx_list,
				connectivity,
				current_state,
				max_coding_degree,
				max_inference_simu,
				effective_max_inflight,
				inactive_nodes,
			)

			for branching_action in branching_actions
				branch_state = deepcopy(current_state)
				push!(branch_state.schedule, branching_action)
				evaluate_state_sim!(connectivity, effective_max_inflight, max_inference_simu, branch_state, solution_knowledge)
				if isinf(branch_state.dist)
					continue
				end

				branch_lower_bound = _astar_inf_v8_scaled_lower_bound(branch_state)
				if branch_lower_bound >= best_scaled_cost
					continue
				end

				branch_key = hash(branch_state)
				if branch_lower_bound >= get(best_lower_by_key, branch_key, Inf)
					continue
				end

				next_state_id += 1
				branch_state.f_score = Float32(branch_lower_bound)
				branch_state.g_score = Float32(branch_lower_bound)
				state_by_id[next_state_id] = branch_state
				state_key_by_id[next_state_id] = branch_key
				parent_by_id[next_state_id] = current_id
				action_by_id[next_state_id] = deepcopy(branching_action)
				best_lower_by_key[branch_key] = branch_lower_bound
				open_set_pq[next_state_id] = branch_lower_bound

				if !isnothing(trace_io) && iteration % trace_every == 0
					_astar_inf_v8_write_trace!(trace_io, "generate", iteration, next_state_id, current_id, branch_state, branch_lower_bound, "", length(open_set_pq), length(closed_lower_by_key), branching_action)
				end
				_astar_inf_v8_maybe_push_snapshot!(snapshots, next_state_id, current_id, branch_state, branch_lower_bound, branching_action)
			end
		end
	finally
		if !isnothing(trace_io)
			close(trace_io)
		end
	end

	snapshot_path = _astar_inf_v8_save_snapshots(snapshots, actual_trace_path)
	if !ismissing(best_solution)
		best_solution.f_score = Float32(best_scaled_cost)
		best_solution.g_score = Float32(best_scaled_cost)
		@info "A* v8 exact traversal finished" best_aoi=_astar_inf_v8_aoi(best_solution) schedule_length=length(best_solution.schedule) expanded=iteration generated=length(state_by_id) trace=actual_trace_path snapshots=snapshot_path
	else
		@warn "A* v8 exact traversal found no terminal state" expanded=iteration generated=length(state_by_id) trace=actual_trace_path snapshots=snapshot_path
	end
	return best_solution
end

function a_star_inf_v8_search(
	connectivity;
	max_coding_degree = 1,
	max_inflight = size(connectivity, 1),
	max_inference_simu = 2,
	solution_knowledge = Dict("t_star" => 0, "dissemination_delays" => zeros(size(connectivity, 1), size(connectivity, 1))),
	star_epsilon = 1.0,
	save_traversal = false,
	trace_path = nothing,
	trace_every = 1,
	save_state_snapshots = false,
	max_expansions = typemax(Int),
)
	return a_star_inf_v8_search(
		connectivity,
		max_coding_degree,
		max_inflight,
		max_inference_simu,
		solution_knowledge,
		star_epsilon;
		save_traversal = save_traversal,
		trace_path = trace_path,
		trace_every = trace_every,
		save_state_snapshots = save_state_snapshots,
		max_expansions = max_expansions,
	)
end
