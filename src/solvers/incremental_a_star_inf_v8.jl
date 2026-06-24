#
# AoI-aligned repair wrapper for incremental_a_star_inf_v3.
#
# The v3 incremental solver keeps the original rolling A* decomposition, but
# its merge score is not fully aligned with the final average-AoI objective.
# Version v8 keeps the v3 search as a candidate and adds a deterministic
# AoI-first repair constructor.  The returned state is the better state under
# the project convention
#
#     average AoI = T / 2 + mean first-packet dissemination delay.
#

using Distributed
include(joinpath(@__DIR__, "incremental_a_star_inf_v3.jl"))

@everywhere begin
	using DrWatson
	quickactivate("age-optimal-multisource-flooding")
	using LinearAlgebra
	using Random

	const _V8_SCORE = (
		age_slope = 0.03,
		launch_factor = 0.20,
		age_receiver = 2.0,
	)

	function _v8_neighbors(connectivity, node_id)
		return findall(connectivity[node_id, :] .> 0)
	end

	function _v8_conflict_matrix(connectivity)
		node_count = size(connectivity, 1)
		conflict = falses(node_count, node_count)
		neighbor_sets = [Set(_v8_neighbors(connectivity, node_id))
						 for node_id in 1:node_count]
		for left in 1:node_count-1
			for right in left+1:node_count
				adjacent = connectivity[left, right] > 0 || connectivity[right, left] > 0
				common_receiver = !isempty(intersect(neighbor_sets[left], neighbor_sets[right]))
				if adjacent || common_receiver
					conflict[left, right] = true
					conflict[right, left] = true
				end
			end
		end
		return conflict
	end

	function _v8_unit_payload(node_count, source_id)
		payload = zeros(Int8, node_count)
		payload[source_id] = 1
		return payload
	end

	function _v8_new_receivers(connectivity, informed, source_id, transmitter)
		return [receiver for receiver in _v8_neighbors(connectivity, transmitter)
				if !informed[source_id, receiver]]
	end

	function _v8_inflight_after_slot(informed, generation_times, selected, time, node_count)
		next_informed = copy(informed)
		next_generation = copy(generation_times)
		for candidate in selected
			source_id = candidate.source
			transmitter = candidate.transmitter
			if next_generation[source_id] == 0
				transmitter == source_id || return typemax(Int)
				next_generation[source_id] = time
			end
		end
		for candidate in selected
			source_id = candidate.source
			for receiver in candidate.receivers
				next_informed[source_id, receiver] = true
			end
		end
		return count(source_id -> next_generation[source_id] != 0 &&
								 !all(next_informed[source_id, :]),
					 1:node_count)
	end

	function _v8_age_first_candidates(connectivity, informed, generation_times, time, rng, noise_scale)
		node_count = size(connectivity, 1)
		candidates = NamedTuple[]
		for transmitter in 1:node_count
			best = nothing
			best_score = -Inf
			for source_id in 1:node_count
				informed[source_id, transmitter] || continue
				receivers = _v8_new_receivers(connectivity, informed, source_id, transmitter)
				isempty(receivers) && continue
				outstanding = node_count - count(informed[source_id, :])
				launched = generation_times[source_id] != 0
				in_flight_age = launched ? time - generation_times[source_id] : 0
				age_pressure = launched ?
					outstanding * (1.0 + _V8_SCORE.age_slope * in_flight_age) :
					_V8_SCORE.launch_factor * outstanding
				score = age_pressure + _V8_SCORE.age_receiver * length(receivers)
				score += noise_scale * randn(rng)
				if score > best_score
					best_score = score
					best = (source = source_id,
							transmitter = transmitter,
							receivers = receivers,
							score = score)
				end
			end
			best === nothing || push!(candidates, best)
		end
		return candidates
	end

	function _v8_select_noninterfering(candidates, conflict, informed,
									  generation_times, time, max_inflight,
									  max_inference_simu)
		node_count = size(conflict, 1)
		ordered = sort(candidates, by = candidate -> candidate.score, rev = true)
		selected = NamedTuple[]
		blocked = falses(node_count)
		for candidate in ordered
			transmitter = candidate.transmitter
			blocked[transmitter] && continue
			tentative = vcat(selected, [candidate])
			_v8_inflight_after_slot(informed, generation_times, tentative, time,
									node_count) > max_inflight && continue
			push!(selected, candidate)
			length(selected) >= max_inference_simu && break
			blocked[transmitter] = true
			for other in 1:node_count
				conflict[transmitter, other] && (blocked[other] = true)
			end
		end
		return selected
	end

	function _v8_age_first_schedule(connectivity; max_inflight = 2,
								   max_inference_simu = 2,
								   starts = 500,
								   seed = 10_000,
								   noise_scale = 0.65,
								   maximum_slots = 1_000_000)
		node_count = size(connectivity, 1)
		conflict = _v8_conflict_matrix(connectivity)
		best = nothing
		for iteration in 1:starts
			rng = MersenneTwister(seed + iteration)
			local_noise = iteration == 1 ? 0.0 : noise_scale
			informed = falses(node_count, node_count)
			generation_times = zeros(Int, node_count)
			arrival_times = zeros(Int, node_count, node_count)
			for source_id in 1:node_count
				informed[source_id, source_id] = true
			end
			schedule = Vector{Vector{Tuple{Int, Int}}}()
			failed = false
			time = 0
			while !all(informed)
				time += 1
				if time > maximum_slots
					failed = true
					break
				end
				candidates = _v8_age_first_candidates(connectivity, informed,
													  generation_times, time,
													  rng, local_noise)
				selected = _v8_select_noninterfering(candidates, conflict, informed,
													 generation_times, time,
													 max_inflight,
													 max_inference_simu)
				if isempty(selected)
					failed = true
					break
				end
				slot = Tuple{Int, Int}[]
				for candidate in selected
					source_id = candidate.source
					transmitter = candidate.transmitter
					if generation_times[source_id] == 0
						transmitter == source_id || error("relay before launch")
						generation_times[source_id] = time
					end
					push!(slot, (transmitter, source_id))
				end
				for candidate in selected
					source_id = candidate.source
					for receiver in candidate.receivers
						if !informed[source_id, receiver]
							informed[source_id, receiver] = true
							arrival_times[source_id, receiver] = time + 1
						end
					end
				end
				push!(schedule, slot)
			end
			failed && continue
			delay_sum = sum(arrival_times[source_id, receiver] -
							generation_times[source_id]
							for source_id in 1:node_count,
								receiver in 1:node_count if source_id != receiver)
			average_aoi = length(schedule) / 2 +
						  delay_sum / (node_count * (node_count - 1))
			result = (schedule = schedule,
					  generation_times = generation_times,
					  arrival_times = arrival_times,
					  delay_sum = delay_sum,
					  average_aoi = average_aoi,
					  iteration = iteration)
			if best === nothing || result.average_aoi < best.average_aoi
				best = result
			end
		end
		best === nothing && error("v8 AoI repair failed to construct a schedule")
		return best
	end

	function _v8_payload_slot_to_action(decoders, payload_slot)
		node_count = size(decoders, 1)
		if length(payload_slot) == 1
			transmitter, source_id = payload_slot[1]
			action = coded_payload_to_action(decoders, transmitter,
											 _v8_unit_payload(node_count, source_id))
			ismissing(action) && error("No action for ($(transmitter),$(source_id))")
			return Dict("at" => "g", "a" => action)
		end
		actions = Tuple{Int, Int}[]
		for (transmitter, source_id) in payload_slot
			action = coded_payload_to_action(decoders, transmitter,
											 _v8_unit_payload(node_count, source_id))
			ismissing(action) && error("No action for ($(transmitter),$(source_id))")
			push!(actions, action)
		end
		return Dict("at" => "m", "a" => Tuple(actions))
	end

	function _v8_replay_payload_schedule(connectivity, payload_schedule, solution_knowledge)
		node_count = size(connectivity, 1)
		start_decoders = zeros(Int8, node_count, node_count, node_count)
		for node_id in 1:node_count
			start_decoders[node_id, 1, node_id] = 1
		end
		state = IncrementalAStarInfState(
			Dict{String, Any}[],
			Dict{String, Any}[],
			start_decoders,
			fill(Int16(32767), node_count, node_count),
			fill(Int16(32767), node_count),
			falses(node_count),
			0.0,
			0.0,
			0.0,
			0.0,
			false,
			Int8[],
			Dict([source_id => source_id for source_id in 1:node_count]),
		)
		for payload_slot in payload_schedule
			action = _v8_payload_slot_to_action(state.decoders, payload_slot)
			push!(state.schedule, action)
			evaluate_incremental_state_inf!(connectivity, state, solution_knowledge)
			state.g_score += state.dist
			state.f_score = state.g_score + state.h
		end
		state.end_state = all(state.content_purged)
		return state
	end

	function _v8_state_aoi(state)
		node_count = length(state.first_encodings)
		first_decodings = copy(state.decoding_levels)
		for node_id in 1:node_count
			first_decodings[node_id, node_id] = 1
		end
		if any(state.first_encodings .== Int16(32767)) ||
		   any(first_decodings .== Int16(32767))
			return Inf
		end
		delay_sum = sum(first_decodings[source_id, receiver] -
						state.first_encodings[source_id]
						for source_id in 1:node_count,
							receiver in 1:node_count if source_id != receiver)
		return length(state.payload_schedule) / 2 +
			   delay_sum / (node_count * (node_count - 1))
	end

	function incremental_a_star_inf_v8_search(connectivity, max_coding_degree,
											 max_inflight, max_inference_simu,
											 solution_knowledge,
											 star_epsilon = 0.0;
											 repair_starts = 500,
											 repair_seed = 10_000,
											 repair_noise_scale = 0.65,
											 maximum_slots = 1_000_000,
											 run_base = true)
		@info "e-inf v8"
		base_state = nothing
		base_aoi = Inf
		if run_base
			base_state = incremental_a_star_inf_search(
				connectivity, max_coding_degree, max_inflight,
				max_inference_simu, deepcopy(solution_knowledge), star_epsilon)
			base_aoi = _v8_state_aoi(base_state)
		end

		repair = _v8_age_first_schedule(
			connectivity;
			max_inflight = max_inflight,
			max_inference_simu = max_inference_simu,
			starts = repair_starts,
			seed = repair_seed,
			noise_scale = repair_noise_scale,
			maximum_slots = maximum_slots,
		)
		repair_state = _v8_replay_payload_schedule(
			connectivity, repair.schedule, deepcopy(solution_knowledge))
		repair_aoi = _v8_state_aoi(repair_state)

		if repair_aoi <= base_aoi
			repair_state.h = Float32(repair_aoi)
			return repair_state
		end
		return base_state
	end
end
