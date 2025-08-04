#
#  Copyright (c) 2023 Institute of Communication Networks (ComNets),
#                     Hamburg University of Technology (TUHH),
#                     https://www.tuhh.de/comnets
#  Copyright (c) 2023 Leonard Fisser <leonard.fisser@tuhh.de>
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 

using LinearAlgebra
using DataStructures
include("$(srcdir())/utilities/bin_mat_rref.jl")
include("$(srcdir())/utilities/aoi_utilities_src.jl")
include("$(srcdir())/utilities/state_utilities_src.jl")
include("$(srcdir())/utilities/topology_utilities.jl")

mutable struct AStarState
	schedule::Vector{Tuple{Int64, Int64}} # The sequence of actions
	payload_schedule::Vector{Tuple{Int64, Vector{Int8}}} # Also sequence of actions, but with binary payload vector
	decoders::Array{Int8, 3} # The decoder matrices of every node
	decoding_levels::Matrix{Int8} # Logs the first time data was decoded at every node -> matrix
	first_encodings::Vector{Int8} # Logs the very first transmissions of data -> vector
	content_purged::BitVector # Data fully disseminated is removed, so we keep track here

	f_score::Float32 # a_star variable
	g_score::Float32# a_star variable
	h::Float32 # a_star variable
	dist::Float32 # a_star variable

	end_state::Bool # Schedule end
	present_content::Vector{Int8} # Data which is currently inflight / is available to be encoded
end

Base.deepcopy(m::AStarState) = AStarState(copy(m.schedule),
	copy(m.payload_schedule),
	copy(m.decoders),
	copy(m.decoding_levels),
	copy(m.first_encodings),
	copy(m.content_purged),
	copy(m.f_score),
	copy(m.g_score),
	copy(m.h),
	copy(m.dist),
	copy(m.end_state),
	copy(m.present_content))

Base.show(io::IO, x::AStarState) = print(io, "F: $(x.f_score), G: $(x.g_score), Schedule: $(x.schedule)")
Base.hash(x::AStarState) = hash(deepcopy(x.decoders))

function evaluate_state!(connectivity, max_inflight, state, solution_knowledge)
	N = size(connectivity, 1)
	step = size(state.schedule, 1)
	decoders = state.decoders
	good_action = false # ACtion changed something and was not a no-op

	if step != 0 #  First state, we don't need to do a full analysis     
		action = state.schedule[step]
		tx_id = action[1]
		good_action, = apply_action!(connectivity, decoders, action)
		if good_action
			payload = action_to_coded_payload(decoders, action)
			push!(state.payload_schedule, (tx_id, payload))
			if state.first_encodings[tx_id] == 127 && payload[tx_id] == 1
				state.first_encodings[tx_id] = step
			end
		end
	end

	dist = 0.0 # Cost associated with this action
	h_dist = 0 #  Remaining expected costs

	is_decoded = get_is_decoded(decoders)
	decoded_indices = findall(is_decoded .== 1)
	for index in decoded_indices
		if state.decoding_levels[index] == 127 # 127==initialisation -> Only log the first decoding event
			state.decoding_levels[index] = step + 1
			dist += 1 # The 1 dissemination delay associated with the decoded symbol
		end
	end

	purge_disseminated_contents!(decoders, state.content_purged) # Remove decoded rows in decoder matrices to simplify rref calculation

	# Actions which resulted in no-op are useless
	if (!good_action) || (length(get_current_inflight_data(step, state.first_encodings, state.decoding_levels)) > max_inflight)
		dist = Inf
	else
		min_t_star = solution_knowledge["t_star"]
		dissemination_delays = solution_knowledge["dissemination_delays"]

		for c_id ∈ findall(.!state.content_purged)
			if state.first_encodings[c_id] < 127
				decoding_status = is_decoded[c_id, :]
				# dist += 0
				dist += (N - sum(decoding_status)) # Inflight, ages additionaly, for every non-delivered-to node
				# We will at least collect the dissemination_delay, for the data which was not intially send yet
				# if (N - sum(decoding_status)) != 0
				#     println("dist add")
				# end

				for d_id ∈ 1:N
					if decoding_status[d_id] == 0
						h_dist += max(0, dissemination_delays[c_id, d_id] - (step - state.first_encodings[c_id]))
					end
				end
			else
				# we did not yet encode, and we will at least collect the content_ids dissemination delay as dist
				h_dist += sum(dissemination_delays[c_id, :]) - dissemination_delays[c_id, c_id]
			end
		end

		if sum(state.decoders) == 0 # all decoded, end state reached
			state.end_state = true
			state.h = 0
		else # estimate how much more aoi cost will be collected until the end-state, should always underestimate
			min_tx_left = max(1, min_t_star - step)
			state.h = (N * (N - 1) / 2) * min_tx_left + h_dist
		end
	end

	dist += (N * (N - 1) / 2) # Taking an action, increases the schedule length by 1 (see age equation for specific factor) # 1/2
	state.dist = dist

	return nothing
end

function a_star_search(connectivity, max_coding_degree, max_inflight, solution_knowledge, star_epsilon = 1.0)
	N = size(connectivity, 1)
	_, non_leaves, true_leaves, pseudo_leaves, mcds = get_node_types(connectivity)

	action_set = vec(collect(Iterators.product(1:N, 1:2^N)))
	start_decoders = zeros(Int8, N, N, N)
	[start_decoders[n, 1, n] = 1 for n ∈ collect(1:N)]
	start = AStarState([], [], (start_decoders), ones(Int8, N, N) .* 127, ones(Int8, N) .* 127, falses(N), 0.0, 0.0, 0.0, 0.0, false, [])
	evaluate_state!(connectivity, max_inflight, start, solution_knowledge)

	step = 1
	open_set_pq = PriorityQueue(hash(start) => start.f_score)
	open_set = Dict(hash(start) => start)

	best_solution = missing

	while !isempty(open_set_pq)

		best = dequeue!(open_set_pq)
		candidate = open_set[best]

		if candidate.end_state
			best_solution = candidate
			break
		end

		# if mod(step, 5000) == 0
		# 	println("Heap size: $(length(open_set_pq)), Dict size: $(length(open_set))")
		# end
		# Every leave node should only transmitted once (its own data), afterwards deactivated
		inactive_nodes = intersect(true_leaves, [action[1] for action in candidate.schedule])
		# All possible actions (includes no-op)
		branching_actions = valid_actions(action_set, candidate, max_coding_degree, inactive_nodes)

		# grow set with all possibilities
		for branching_action in branching_actions
			branch_candidate = deepcopy(candidate)
			branch_candidate.g_score = Inf

			# Add action to schedule
			push!(branch_candidate.schedule, branching_action)
			evaluate_state!(connectivity, max_inflight, branch_candidate, solution_knowledge)

			if !isinf(branch_candidate.dist)
				tentative_g_score = candidate.g_score + branch_candidate.dist
				check = get(open_set, hash(branch_candidate), missing)

				if ismissing(check)
					check = branch_candidate
				end

				if tentative_g_score < check.g_score
					check = deepcopy(branch_candidate)
					check.g_score = tentative_g_score
					check.f_score = tentative_g_score + star_epsilon * check.h
					open_set[hash(check)] = check
					open_set_pq[hash(check)] = check.f_score
				end
			end
		end
		step += 1
	end
	# GC.gc()  # 手动触发垃圾回收
	return best_solution
end
