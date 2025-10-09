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
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 这里把struct Base 和 估计状态以及估计合并作为 设置为everywhere
@everywhere begin
	using DrWatson
	quickactivate("age-optimal-multisource-flooding")
	using DataStructures
	include("$(srcdir())/utilities/aoi_utilities_src_inf_v3.jl")
	mutable struct IncrementalAStarInfState
		schedule::Vector{Dict{String, Any}} # The sequence of actions
		payload_schedule::Vector{Dict{String, Any}} # Also sequence of actions, but with binary payload vector
		decoders::Array{Int8, 3} # The decoder matrices of every node
		decoding_levels::Matrix{Int16} # Logs the first time data was decoded at every node -> matrix
		first_encodings::Vector{Int16} # Logs the very first transmissions of data -> vector
		content_purged::BitVector # Data fully disseminated is removed, so we keep track here
		f_score::Float32 # a_star variable
		g_score::Float32 # a_star variable
		h::Float32 # a_star variable
		dist::Float32 # a_star variable
		end_state::Bool # Schedule end
		present_content::Vector{Int8}
		active_content_map::Dict{Int64, Int64} # 添加了一个映射的map
		# max_inference::Int8
	end
	Base.deepcopy(m::IncrementalAStarInfState) = IncrementalAStarInfState(copy(m.schedule),
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
		copy(m.present_content),
		copy(m.active_content_map))
	# copy(m.max_inference))
	Base.show(io::IO, x::IncrementalAStarInfState) = print(io, "F: $(x.f_score), G: $(x.g_score), Schedule: $(x.schedule)")
	Base.hash(x::IncrementalAStarInfState) = hash(deepcopy(x.decoders))
	function evaluate_incremental_state_inf!(connectivity, state, solution_knowledge)
		N = size(connectivity, 1)
		step = size(state.schedule, 1)
		decoders = state.decoders
		good_action_all_simu = false
		if step != 0
			action_s_simu = state.schedule[step]
			good_action_all_simu, UpdateNeibors_all_simu_dict = apply_action_inf!(connectivity, decoders, action_s_simu)
			if good_action_all_simu
				reduced_payload_s_simu = action_to_coded_payload_all_simu(decoders, action_s_simu)
				if action_s_simu["at"] == "g"
					action = action_s_simu["a"]
					tx_id = action[1]
					reduced_payload = reduced_payload_s_simu["p"]
					payload = zeros(Int8, N)
					for x in state.active_content_map
						if x[2] != 0
							@inbounds if reduced_payload[x[1]] == 1
								@inbounds payload[x[2]] = 1
							end
						end
					end
					payload_dict_after_judge_reduced = Dict(
						"pt" => "g",
						"t" => action[1],
						"p" => payload,
					)
					push!(state.payload_schedule, payload_dict_after_judge_reduced) # 
					@inbounds if (state.first_encodings[tx_id] == typemax(typeof(state.first_encodings[tx_id]))) && (payload[tx_id] == 1)
						@inbounds state.first_encodings[tx_id] = step
					end
				elseif action_s_simu["at"] == "m"
					payload_s_dict_after_judge_reduced = []
					payload_s_simu_temp = reduced_payload_s_simu["p"]
					tx_list = reduced_payload_s_simu["t"]
					action_s = action_s_simu["a"]
					index_1 = 1 # 第几个action，对应第几个payload
					for action in action_s # 同时传输的动作
						tx_id = action[1]
						reduced_payload = payload_s_simu_temp[index_1] #同时生成的负载
						payload = zeros(Int8, N)
						for x in state.active_content_map
							if x[2] != 0
								@inbounds if reduced_payload[x[1]] == 1
									@inbounds payload[x[2]] = 1
								end
							end
						end
						push!(payload_s_dict_after_judge_reduced, payload)
						@inbounds if (state.first_encodings[tx_id] == typemax(typeof(state.first_encodings[tx_id]))) && (payload[tx_id] == 1)
							@inbounds state.first_encodings[tx_id] = step
						end
						index_1 = index_1 + 1
					end
					payload_dict_after_judge_reduced = Dict(
						"pt" => "m",
						"t" => tx_list,
						"p" => payload_s_dict_after_judge_reduced,
					)
					push!(state.payload_schedule, payload_dict_after_judge_reduced) # 
				end
			else
				state.dist = Inf
				return nothing
			end
		end
		dist = 0.0 # Cost associated with this action
		h_dist = 0 # Remaining expected costs
		is_decoded = get_is_decoded(decoders, state.active_content_map)# 这个函数出现了重新构造
		decoded_indices = findall(is_decoded .== 1) # 找出所有已解码的索引
		for index in decoded_indices
			@inbounds if state.decoding_levels[index[1], index[2]] == 32767
				#TODO: improve magic number 改进 magic number 
				# 程式设计中所谓的直接写在程式码里的具体数值
				@inbounds state.decoding_levels[index[1], index[2]] = step + 1 #记录解码时间
				dist += 1
			end
		end

		before_purge = deepcopy(state.content_purged)
		purge_disseminated_contents!(decoders, state.content_purged, state.active_content_map)
		if (sum(before_purge) < sum(state.content_purged))
			a_p = findall(state.content_purged .!= before_purge)
			for a in a_p
				for x in state.active_content_map
					@inbounds if x[2] == a
						@inbounds state.active_content_map[x[1]] = 0
					end
				end
			end
		end

		# dissemination_delays = solution_knowledge["dissemination_delays"]
		# for c_id ∈ findall(.!state.content_purged)
		# 	@inbounds if state.first_encodings[c_id] < 32767 #TODO:
		# 		@inbounds decoding_status = is_decoded[c_id, :]
		# 		# dist += 0
		# 		# dist += (N - sum(decoding_status))  # 只在nc(network coding) 方案中用, sys,sys-e,sys-inf,sys-inf-e都用不到
		# 		# Inflight, ages additionaly, for every non-delivered-to node
		# 		# We will at least collect the dissemination_delay, so we should consider this in the heuristic
		# 		# 对于每个未交付的节点，在飞行中，年龄增加
		# 		# 我们至少会收集dissemination_delay，用于尚未发送的数据
		# 		# for d_id ∈ 1:N
		# 		# 	@inbounds if decoding_status[d_id] == 0
		# 		# 		@inbounds h_dist += max(0, dissemination_delays[c_id, d_id] - (step - state.first_encodings[c_id]))
		# 		# 	end
		# 		# end
		# 	else
		# 		# we did not yet encode, and we will at least collect the content_ids dissemination delay as dist
		# 		# 我们还没有编码，我们至少会收集传播延迟content_ids作为 dist
		# 		# @inbounds h_dist += sum(dissemination_delays[c_id, :]) - dissemination_delays[c_id, c_id]
		# 	end
		# end

		if sum(state.decoders) == 0
			state.end_state = true
			state.h = 0
		else
			# min_tx_left = max(1, solution_knowledge["t_star"] - step)
			# state.h = (N * (N - 1) / 2) * min_tx_left + h_dist
			state.h = 0
		end
		# Add additional aging due to increase in update interval
		# dist += (N * (N - 1) / 2) # T带来的
		state.dist = dist
		return nothing
	end


	# 评估合并
	function evaluate_merge_inf(loop_state, merge, action_set, simu_tx_list, connectivity, max_coding_degree, max_inflight, max_inference_simu, solution_knowledge, t_star, star_epsilon, true_leaves)
		N = size(connectivity, 1)
		reduced_loop_state = deepcopy(loop_state)
		reduced_loop_state.decoders = zeros(Int8, max_inflight, max_inflight, N)
		reduced_loop_state.end_state = false
		solution_knowledge["t_star"] = t_star
		new_content_map = Dict([x => 0 for x in 1:max_inflight])
		union!(reduced_loop_state.present_content, merge)
		for content_id in reduced_loop_state.present_content
			@inbounds if sum(loop_state.decoders[content_id, :, :]) > 0
				# find unused dict entry
				new_entry = missing
				for x in new_content_map
					@inbounds if x[2] == 0
						@inbounds new_content_map[x[1]] = content_id
						new_entry = x[1]
						break
					end
				end
				node_ranks = get_node_ranks(reduced_loop_state.decoders)
				for node_id ∈ 1:N
					@inbounds row_entry = node_ranks[node_id] + 1
					for row_id ∈ 1:N
						@inbounds if loop_state.decoders[content_id, row_id, node_id] == 1
							@inbounds reduced_loop_state.decoders[new_entry, row_entry, node_id] = loop_state.decoders[content_id, row_id, node_id]
							row_entry += 1
						end
					end
				end
			end
		end
		reduced_loop_state.active_content_map = deepcopy(new_content_map)
		step = 1
		open_set_pq = PriorityQueue(hash(reduced_loop_state) => reduced_loop_state.f_score)
		open_set = Dict(hash(reduced_loop_state) => reduced_loop_state)
		best_solution = missing
		while !isempty(open_set_pq)
			best = dequeue!(open_set_pq)
			candidate = open_set[best]
			if candidate.end_state
				best_solution = deepcopy(candidate)
				break
			end
			branching_actions = valid_actions_all_simu(action_set, simu_tx_list, connectivity, candidate, max_coding_degree, max_inference_simu, max_inflight, [])
			for branching_action in branching_actions
				branch_candidate = deepcopy(candidate)
				branch_candidate.g_score = Inf
				push!(branch_candidate.schedule, branching_action)
				evaluate_incremental_state_inf!(connectivity, branch_candidate, solution_knowledge)
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
		return best_solution
	end
end

########### 这两个函数不参与并行计算
function get_content_costs(connectivity, simu_tx_list)
	N = size(connectivity, 1)
	solution_knowledge = Dict([("dissemination_delays" => zeros(N, N)), ("t_star", 0)])
	content_costs = zeros(N)
	action_set = vec(collect(Iterators.product(1:N, 1)))
	sys_merges = collect(enumerate(collect(1:N)))
	start_decoders = zeros(Int8, N, N, N)
	[start_decoders[n, 1, n] = 1 for n ∈ collect(1:N)]
	sys_state = IncrementalAStarInfState([], [], start_decoders, ones(Int16, N, N) .* 32767, ones(Int16, N) .* 32767, falses(N), 0.0, 0.0, 0.0, 0.0, false, [], Dict([x => 0 for x in 1:1]))
	sys_combination_log = map(merge -> evaluate_merge_inf(deepcopy(sys_state), merge[1], action_set, simu_tx_list, connectivity, 1, 1, 1, solution_knowledge, 0, 0.0, []), sys_merges)
	for x in sys_combination_log
		content_costs[x.present_content[1]] = x.g_score
	end
	return content_costs
end
function incremental_a_star_inf_search(connectivity, max_coding_degree, max_inflight, max_inference_simu, solution_knowledge, star_epsilon = 0.0)
	@info "e-inf v3"
	N = size(connectivity, 1)
	start_decoders = zeros(Int8, N, N, N)
	[start_decoders[n, 1, n] = 1 for n ∈ collect(1:N)]
	accumulated_state = IncrementalAStarInfState([], [], start_decoders, ones(Int16, N, N) .* 32767, ones(Int16, N) .* 32767, falses(N), 0.0, 0.0, 0.0, 0.0, false, [], Dict([x => 0 for x in 1:max_inflight]))
	simu_tx_list = []
	if max_inference_simu > 1
		N_t = UInt32(N)
		all_C = collect(1:(2^N_t-1))
		for simu_tx_node_set in all_C
			node_list = get_all_index(simu_tx_node_set)
			n = length(node_list) # 这里面是一组节点
			flage = false
			for i in 1:n-1
				@inbounds i_neighbor = findall(connectivity[node_list[i], :] .> 0)
				for j in i+1:n
					@inbounds j_neighbor = findall(connectivity[node_list[j], :] .> 0)
					@inbounds if length(intersect(j_neighbor, i_neighbor)) == 0 && length(intersect(j_neighbor, [node_list[i]])) == 0 && length(intersect(i_neighbor, [node_list[j]])) == 0
						flage = true
					else
						flage = false
						break
					end
				end
				if !flage
					break
				end
			end
			if flage
				if length(node_list) <= max_inference_simu
					push!(simu_tx_list, node_list)
				end
			end
		end
	end
	content_costs = get_content_costs(connectivity, simu_tx_list)
	action_set = vec(collect(Iterators.product(1:N, 1:(2^max_inflight)-1)))
	true_leaves = []
	cached_merges = []
	minimum_purges = 1
	while sum(accumulated_state.content_purged) < N
		# Calculate all possibile content merges
		num_draws = min(max_inflight - (length(accumulated_state.present_content) - sum(accumulated_state.content_purged)), N - length(accumulated_state.present_content))
		draw_from = setdiff(collect(1:N), accumulated_state.present_content)
		merges = []
		if (num_draws == max_inflight) && !isempty(cached_merges)
			for merge in cached_merges
				@inbounds if !any([m in accumulated_state.present_content for m in merge[2]])
					@inbounds merges = [(1, deepcopy(merge[2]))]
					break
				end
			end
		end
		if isempty(merges)
			merges = collect(enumerate(combinations(draw_from, num_draws)))
		end
		t_stars = zeros(size(merges, 1))
		combination_log =
			pmap(merge -> evaluate_merge_inf(deepcopy(accumulated_state), merge[2], action_set, simu_tx_list, connectivity, max_coding_degree, max_inflight, max_inference_simu, solution_knowledge, t_stars[merge[1]], star_epsilon, true_leaves), merges)
		sorter = sortperm(combination_log, by = v -> sum(content_costs[v.present_content]) - v.g_score, rev = true)
		@inbounds combination_log = combination_log[sorter]
		merges = merges[sorter]
		if isempty(cached_merges)
			cached_merges = deepcopy(merges)
		end
		best_partial_state = deepcopy(combination_log[1])
		roll_back_state = IncrementalAStarInfState([], [], zeros(Int8, N, N, N), ones(Int16, N, N) .* 32767, ones(Int16, N) .* 32767, falses(N), 0.0, 0.0, 0.0, 0.0, false, [], Dict([x => x for x in 1:N]))
		roll_back_state.present_content = deepcopy(best_partial_state.present_content)
		for x ∈ 1:N
			@inbounds roll_back_state.decoders[x, 1, x] = 1
		end
		evaluate_incremental_state_inf!(connectivity, roll_back_state, solution_knowledge)
		for pay_sim in best_partial_state.payload_schedule  # 这个是payload 
			if sum(roll_back_state.content_purged) >= minimum_purges
				break
			end
			# payload_simu = deepcopy(pay_sim)
			if pay_sim["pt"] == "g"
				tx_id = pay_sim["t"]
				coded_payload = pay_sim["p"] # 负载
				# 负载到动作
				roll_back_action = coded_payload_to_action(roll_back_state.decoders, tx_id, coded_payload)
				roll_back_action_dict = Dict(
					"at" => "g",
					"a" => roll_back_action
				)
				# push!(roll_back_state.schedule, roll_back_action)
				push!(roll_back_state.schedule, roll_back_action_dict)

			elseif pay_sim["pt"] == "m"
				sim_action = []
				payload_s = pay_sim["p"]
				tx_list = pay_sim["t"]
				index_2 = 1 # 第几个action，对应第几个payload
				for coded_payload in payload_s # 负载
					@inbounds tx_id = tx_list[index_2]
					roll_back_action = coded_payload_to_action(roll_back_state.decoders, tx_id, coded_payload)
					push!(sim_action, roll_back_action)
					index_2 = index_2 + 1
				end
				roll_back_action_dict = Dict(
					"at" => "m",
					"a" => sim_action,
					# "sim_tx_number" => length(sim_action)
				)
				# push!(roll_back_state.schedule, roll_back_action)
				push!(roll_back_state.schedule, roll_back_action_dict)
			end
			evaluate_incremental_state_inf!(connectivity, roll_back_state, solution_knowledge)
			tentative_g_score = roll_back_state.g_score + roll_back_state.dist
			roll_back_state.g_score = tentative_g_score
			roll_back_state.f_score = tentative_g_score + star_epsilon * roll_back_state.h
		end
		accumulated_state = deepcopy(roll_back_state)
		minimum_purges = sum(roll_back_state.content_purged) + 1
		# sleep(1)
	end
	accumulated_state.schedule = payload_schedule_to_action_schedule_inf(connectivity, accumulated_state.payload_schedule)
	# GC.gc()  # 手动触发垃圾回收
	return accumulated_state
end
