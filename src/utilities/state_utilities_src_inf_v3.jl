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

include("$(srcdir())/utilities/bin_mat_rref.jl")

function get_all_index(bibdata)
	Index_set = []
	data = bibdata
	get_index = 1
	while data > 0
		temp = data & 1
		if temp == 1
			push!(Index_set, get_index)
		end
		data = data >> 1
		get_index += 1
	end

	return Index_set
end

# function Get_Simu_tx_actions(connectivity)
# 	N = size(connectivity, 1)
# 	N = UInt32(N)
# 	# @info N
# 	all_C = collect(1:(2^N-1))
# 	simu_tx_node_set_list = []
# 	for simu_tx_node_set in all_C
# 		# print("十进制表示的数字为 $simu_tx_node_set ")
# 		node_list = get_all_index(simu_tx_node_set)
# 		n = length(node_list) # 这里面是一组节点
# 		# println("节点集合为 $node_list")
# 		flage = false
# 		for i in 1:n-1
# 			i_neighbor = findall(connectivity[node_list[i], :] .> 0)
# 			for j in i+1:n
# 				j_neighbor = findall(connectivity[node_list[j], :] .> 0)
# 				@inbounds if length(intersect(j_neighbor, i_neighbor)) == 0 && length(intersect(j_neighbor, [node_list[i]])) == 0 && length(intersect(i_neighbor, [node_list[j]])) == 0
# 					flage = true
# 				else
# 					flage = false
# 					break
# 				end
# 			end
# 			if !flage
# 				break
# 			end
# 		end

# 		if flage
# 			push!(simu_tx_node_set_list, node_list)
# 		end
# 		# 然后是让 节点发送的时候不接收
# 	end

# 	# sort!(simu_tx_node_set_list)
# 	return simu_tx_node_set_list
# end

function is_row_encoded(coding_action, row_id)
	"""Returns true, if row_id is part used in the coding_action."""
	return (coding_action >> (row_id - 1) & 1) == 1
end

function get_node_ranks(state)
	N = size(state, 3)
	node_ranks = zeros(Int, N)
	for node in 1:N
		@inbounds node_ranks[node] = sum(sum(state[:, :, node], dims = 1) .!= 0)
	end
	return node_ranks
end


#即便在inf情景中, 同时传输的节点, 每次也只能传输一个信息,并且没有相同的接收
# 因此 node_id 不会是同一个
# 
# 1---> A(2,3,4)---B(1,2,3)<---4
# 信息 1 2 3 4 完成了在所有点上的传输，由于接收节点是不同的
# 因此 sum(state(:row_id,node_id))上的没有影响，node_id不影响
function get_is_decoded(state)
	N = size(state, 3)
	is_decoded = zeros(Int8, N, N)
	for node_id ∈ 1:N
		for row_id ∈ 1:size(state, 2)
			row_sum = sum(state[:, row_id, node_id])
			# if row_sum > 0
			# 	is_recived = findall(state[:, row_id, node_id] .> 0)   # 找到有内容的号
			# 	for index in is_recived
			# 		if index > 0
			# 			is_decoded[index, node_id] = 1
			# 		end
			# 	end

			if row_sum == 1
				@inbounds content_id = findfirst(state[:, row_id, node_id] .== 1)
				@inbounds is_decoded[content_id, node_id] = 1

			elseif row_sum == 0
				break
			end
		end
	end
	return is_decoded
end

function get_is_decoded(state, content_map)
	N = size(state, 3)
	is_decoded = zeros(Int8, N, N)
	for node_id ∈ 1:N
		for row_id ∈ 1:size(state, 2)
			row_sum = sum(state[:, row_id, node_id])
			# if row_sum > 0
			# 	is_recived = findall(state[:, row_id, node_id] .> 0)
			# 	for index in is_recived
			# 		if index > 0
			# 			is_decoded[content_map[index], node_id] = 1
			# 		end
			# 	end
			if row_sum == 1
				@inbounds content_id = findfirst(state[:, row_id, node_id] .== 1)
				@inbounds is_decoded[content_map[content_id], node_id] = 1
			elseif row_sum == 0
				break
			end
		end
	end
	return is_decoded
end

function action_to_coded_payload(state, action)
	N = size(state, 3)
	payload = zeros(Int8, size(state, 1))
	# println("action $action")
	tx_id = action[1]
	coding_action = action[2]
	for coding_row ∈ 1:N
		if is_row_encoded(coding_action, coding_row)
			@inbounds payload .+= state[:, coding_row, tx_id]
		end
	end
	return payload
end

function action_to_coded_payload_inf(state, tx_id, coding_action)
	N = size(state, 3)
	payload = zeros(Int8, size(state, 1))
	# tx_id = action[1]
	# coding_action = action[2]
	for coding_row ∈ 1:N
		if is_row_encoded(coding_action, coding_row)
			@inbounds payload .+= state[:, coding_row, tx_id]
		end
	end
	return payload
end


function action_to_coded_payload_all_simu(state, action_s_simu) # action_s_simu 中是同时传输的动作
	N = size(state, 3)
	# println("state=> $state")
	payload = zeros(Int8, size(state, 1))
	payload_sim = []
	if action_s_simu["at"] == "g"
		action = action_s_simu["a"]
		# println("action $action")
		payload_temp = action_to_coded_payload(state, action)
		payload_sim = Dict(
			"pt" => "g",
			"t" => action[1],
			"p" => payload_temp,
		)
		# println(" action=> $action -> $payload")
		# payload_sim = deepcopy(payload_dict)
		# @info payload_sim

	elseif action_s_simu["at"] == "m"
		action = action_s_simu["a"]
		payload_s = []
		tx_list = []
		# print("action_s $action_s =>")
		for a in action
			# println("action_s $action_s")
			payload_temp = action_to_coded_payload(state, a)
			push!(payload_s, payload_temp)
			push!(tx_list, a[1])
			# print(" $action -> $payload")
			# tx_list =[]
		end
		# println()
		payload_sim = Dict(
			"pt" => "m",
			"t" => tx_list,
			"p" => payload_s
		)
		# payload_sim = deepcopy(payload_dict)
		# @info payload_sim
	end
	# Dict("payload_type"=> "single, "payload" = [1,0,1])
	# Dict("payload_type"=> "sim, "payload_sim" = [[1,0,1],[1,0,0]])
	# println("payload_sim=>$payload_sim")
	return payload_sim
end

function apply_action_inf!(connectivity, state, action_s_simu)
	"""Transforms the input state according to action tuple (tx_id,coding_action) and connectivity."""
	# 进来的动作要进行评估
	# [ Info: Dict{String, Any}("action" => (1, 2), "action_type" => "single")
	# [ Info: Dict{String, Any}("action" => (1, 2), "action_type" => "single")
	# [ Info: Dict{String, Any}("action" => ((1, 2), (3, 4)), "action_type" => "sim")
	N = size(state, 1)
	good_transformation_all_simu = false
	UpdateNeibors_all_simu = []
	good_transformation = false
	UpdateNeibors = []
	connections = []
	payload = []
	UpdateNeibors_all_simu_dict = []
	flage_all_is_good = []

	if action_s_simu["at"] == "g"
		action = action_s_simu["a"]
		tx_id = action[1]
		connections = findall(connectivity[tx_id, :] .> 0)
		# UpdateNeibors = connections
		# println("action $action")
		payload = action_to_coded_payload(state, action)
		for rx_id in connections
			@inbounds zero_row = findfirst(sum(state[:, :, rx_id], dims = 1) .== 0)
			if !isnothing(zero_row)
				zero_row = zero_row[2]
				for c_d ∈ 1:N
					@inbounds state[c_d, zero_row, rx_id] = (state[c_d, zero_row, rx_id] + payload[c_d]) % 2
					#即便是inf情景中,因为没有共同的接收节点，有效的接收中，接受节点每次依然每次接受1个
				end
				bin_mat_rref_line!(state, rx_id, zero_row)
				@inbounds if sum(state[:, zero_row, rx_id]) > 0
					good_transformation = true
				end
			end
		end
		good_transformation_all_simu = good_transformation
		# if good_transformation_all_simu
		# 	UpdateNeibors_all_simu_dict_temp = Dict(
		# 		"neighbor_type" => "single",
		# 		"UpdateNeibors" => UpdateNeibors,
		# 	)
		# 	UpdateNeibors_all_simu_dict = deepcopy(UpdateNeibors_all_simu_dict_temp)
		# end
	elseif action_s_simu["at"] == "m"
		action = action_s_simu["a"] # 一系列动作
		for a in action
			tx_id = a[1]
			good_transformation = false

			@inbounds connections = findall(connectivity[tx_id, :] .> 0)
			# UpdateNeibors = connections
			# println("action_s $action_s")
			payload = action_to_coded_payload(state, a)
			for rx_id in connections
				@inbounds zero_row = findfirst(sum(state[:, :, rx_id], dims = 1) .== 0)
				if !isnothing(zero_row)
					@inbounds zero_row = zero_row[2]
					for c_d ∈ 1:N
						@inbounds state[c_d, zero_row, rx_id] = (state[c_d, zero_row, rx_id] + payload[c_d]) % 2
					end

					bin_mat_rref_line!(state, rx_id, zero_row)
					@inbounds if sum(state[:, zero_row, rx_id]) > 0
						good_transformation = true
						push!(flage_all_is_good, 1)  # 有时候是1个有时候是两个
					end
				end
			end
			if !good_transformation
				# UpdateNeibors_all_simu = []
				good_transformation_all_simu = false
				break
			elseif good_transformation
				# push!(UpdateNeibors_all_simu, UpdateNeibors)
			end
		end
		if length(flage_all_is_good) > 0  # 至少是1个good 邻居多了呢是多个good但是要比 循环数多
			if sum(flage_all_is_good) >= length(action)

				# UpdateNeibors_all_simu_dict_temp = Dict(
				# 	"neighbor_type" => "sim",
				# 	"UpdateNeibors" => UpdateNeibors_all_simu
				# )
				# UpdateNeibors_all_simu_dict = deepcopy(UpdateNeibors_all_simu_dict_temp)
				good_transformation_all_simu = true
			else
				good_transformation_all_simu = false
			end
		else
			good_transformation_all_simu = false
		end
	end
	return good_transformation_all_simu, UpdateNeibors_all_simu_dict
end

function valid_actions(all_actions, state, max_coding_degree, inactive_nodes = [])
	node_ranks = get_node_ranks(state.decoders)
	# 节点发送的信息量<= 节点上已经有的信息量    发送谁的信息<= 飞行中的信息量
	# 以2填充的位数,
	@inbounds step_actions = [action for action in all_actions if (ndigits(action[2], base = 2) <= node_ranks[action[1]]) && (count_ones(action[2]) <= max_coding_degree)]
	# Remove actions of inactive node
	for inactive_node in inactive_nodes
		@inbounds step_actions = [action for action in step_actions if action[1] != inactive_node]
	end
	return step_actions
end


function valid_actions_all_simu(all_actions,simu_tx_list, connectivity, state, max_coding_degree, max_inference_simu, max_inflight, inactive_nodes = []) #v2
	# @info "基本动作集 $all_actions"
	# action_set = vec(collect(Iterators.product(1:N, 1:2^N))) # 基本动作集
	node_ranks = get_node_ranks(state.decoders)
	# @info max_coding_degree

	@inbounds step_actions = [action for action in all_actions if (ndigits(action[2], base = 2) <= node_ranks[action[1]]) && (count_ones(action[2]) <= max_coding_degree)]
	#  step_actions = [action for action in all_actions if (  ndigits(action[2], base = 2) <= node_ranks[action[1]] && ndigits(action[2], base = 2)>0) && (count_ones(action[2]) <= max_coding_degree)] 是不会等于0的

	# 对筛选后的基本动作集, 按照叶子只发一次的规则,去掉已经发过的叶子,Remove actions of inactive node
	for inactive_node in inactive_nodes
		@inbounds step_actions = [action for action in step_actions if action[1] != inactive_node]
	end

	step_actions_single_and_sim = []

	if max_inference_simu > 1

		TxSimuActionsPool = []
		simu_tx_list_new = []

		@inbounds all_can_tx = [i[1] for i in step_actions]
		for tx_s in simu_tx_list
				flage_all_is_in_actin_set = []
				for tx in tx_s
					if length(findall(all_can_tx .== tx)) > 0
						push!(flage_all_is_in_actin_set, 1)
					else
						flage_all_is_in_actin_set = [0]
					end
				end

				if sum(flage_all_is_in_actin_set) == length(tx_s)
					push!(simu_tx_list_new, tx_s)
				end
			# end
		end

		step_actions_sim = []

		for tx_list in simu_tx_list_new # [1,2]
			LIST = []
			for tx in tx_list
				@inbounds temp = findall(all_can_tx[:] .== tx)
				LIST_temp = []

				len_temp = length(temp)
				if (len_temp > 0)
					for i in 1:len_temp
						@inbounds push!(LIST_temp, step_actions[temp[i]]) # 节点可以搞的所有的动作比如 下面的例子 
					end
					temp = [Dict(
						"t" => tx,
						"a" => LIST_temp,
					)]
					union!(TxSimuActionsPool, temp)
				end
			end
		end

		for tx_list in simu_tx_list_new # [1,2]
			LIST = []
			for i in tx_list  
				list1 = []
				flage = false
				for pool in TxSimuActionsPool
					if pool["t"] == i
						list1 = pool["a"]
						flage = true
					end
				end
				if flage
					push!(LIST, list1)
				end
			end

			if length(tx_list) == 2
				@inbounds action_set = vec(collect(Iterators.product(LIST[1], LIST[2])))
				for i in action_set
					push!(step_actions_sim, i)
				end
			elseif length(tx_list) == 3
				@inbounds action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3])))
				for i in action_set
					push!(step_actions_sim, i)
				end
				# elseif max_color_num == 4
				# 	 action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3], LIST[4])))
				# 	# @info "这些动作n*m*m*m' 组合 $action_set  "
				# 	for i in action_set
				# 		push!(step_actions_sim, i)
				# 	end
				# elseif max_color_num == 5
				# 	 action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3], LIST[4], LIST[5])))
				# 	# @info "这些动作n*m*m*m*m组合 $action_set  "
				# 	for i in action_set
				# 		push!(step_actions_sim, i)
				# 	end
				# elseif max_color_num == 6
				# 	 action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3], LIST[4], LIST[5], LIST[6])))
				# 	# @info "这些动作n*m*m*m*m*m组合 $action_set  "
				# 	for i in action_set
				# 		push!(step_actions_sim, i)
				# 	end
				# elseif max_color_num == 7
				# 	 action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3], LIST[4], LIST[5], LIST[6], LIST[7])))
				# 	# @info "这些动作n*m*m*m*m*m组合 $action_set  "
				# 	for i in action_set
				# 		push!(step_actions_sim, i)
				# 	end

				# elseif max_color_num == 8
				# 	action_set = vec(collect(Iterators.product(LIST[1], LIST[2], LIST[3], LIST[4], LIST[5], LIST[6], LIST[7],LIST[8])))
				# 	# @info "这些动作n*m*m*m*m*m组合 $action_set  "
				# 	for i in action_set
				# 		push!(step_actions_sim, i)
				# 	end
			end
			#...... 可以继续写同时多少个并发的
			# else
			# println("not satisify the sim_tx_number band")
			# end

		end

		for sim_action in step_actions_sim
			# 在max_degree=1的情景下不用判断
			# @info "同时传输的动作的信息量的和也不大于max_inflight 主要应对每次传输的信息量不为1的情况$sim_action" 
			degree_tx_sum = 0
			for action in sim_action
				@inbounds degree_tx_sum = degree_tx_sum + count_ones(action[2])
			end
			# @info "degree_tx_sum = $degree_tx_sum"
			if degree_tx_sum <= max_inflight
				sim_type_and_action = Dict(
					"at" => "m",
					"a" => sim_action
					# "sim_tx_number" => length(sim_action)
				)
				# println(sim_type_and_action)
				push!(step_actions_single_and_sim, sim_type_and_action)
			end
		end
	end

	# 单点不并发
	for single_action in step_actions
		single_type_and_action = Dict(
			"at" => "g",
			"a" => single_action
		)
		push!(step_actions_single_and_sim, single_type_and_action)
		# println(single_type_and_action)
	end

	return step_actions_single_and_sim
end


function purge_disseminated_contents!(state)
	content_purged = falses(size(state, 3))
	purge_disseminated_contents!(state, content_purged)
	return nothing
end
function purge_disseminated_contents!(state, content_purged)
	content_map = Dict([x => x for x in 1:size(state, 3)])
	purge_disseminated_contents!(state, content_purged, content_map)
	return nothing
end

function purge_disseminated_contents!(state, content_purged, content_map)
	N = size(state, 3)
	N1 = size(state, 1)
	N2 = size(state, 2)
	is_decoded = get_is_decoded(state)
	for content_id in findall(sum(is_decoded, dims = 2) .== N)
		content_id = content_id[1]
		# if content_id > 0
		content_purged[content_map[content_id]] = true
		for node_id ∈ 1:N
			for row_id ∈ 1:N2
				@inbounds if state[content_id, row_id, node_id] == 1
					@inbounds state[:, row_id, node_id] .= 0
					bin_mat_rref!(state, node_id)
					break
				end
			end
			# end
		end
	end
	return nothing
end


function coded_payload_to_action(state, tx_id, coded_payload)
	N = size(state, 3)
	possible_actions = 0:sum(binomial.(N, 1:N))
	for action in possible_actions
		test_payload = zeros(N)
		for row_id in 1:N
			if is_row_encoded(action, row_id)
				test_payload = (test_payload + state[:, row_id, tx_id])
			end
			if test_payload == coded_payload
				# return (tx_id, action)
				return (tx_id, action)
			end
		end
	end
	return missing
end


# function payload_schedule_to_action_schedule(connectivity, payload_schedule)
# 	N = size(connectivity, 1)
# 	state = zeros(Int8, N, N, N)
# 	for node_id ∈ 1:N
# 		state[node_id, 1, node_id] = 1
# 	end
# 	action_schedule = []
# 	for payload in payload_schedule
# 		# action = coded_payload_to_action(state, payload[1], payload[2])
# 		action = coded_payload_to_action(state, tx_id, payload)  
# 		apply_action_inf!(connectivity, state, action)
# 		purge_disseminated_contents!(state, get_is_decoded(state))
# 		push!(action_schedule, action)
# 	end
# 	return action_schedule
# end

function apply_action!(connectivity, state, action)
	"""Transforms the input state according to action tuple (tx_id,coding_action) and connectivity."""
	N = size(state, 1)
	tx_id = action[1]

	good_transformation = false
	connections = findall(connectivity[tx_id, :] .> 0)
	# UpdateNeibors = connections
	UpdateNeibors = []
	payload = action_to_coded_payload(state, action)
	# println("rx_id judge_________________________________________________________________________")
	# good_num = 0
	for rx_id in connections
		zero_row = findfirst(sum(state[:, :, rx_id], dims = 1) .== 0)
		if !isnothing(zero_row)
			# println("in judge")
			zero_row = zero_row[2]
			for c_d ∈ 1:N
				@inbounds state[c_d, zero_row, rx_id] = (state[c_d, zero_row, rx_id] + payload[c_d]) % 2
			end
			bin_mat_rref_line!(state, rx_id, zero_row)
			if sum(state[:, zero_row, rx_id]) > 0
				good_transformation = true
				# good_num = good_num + 1
				# println("in judge $good_num")
			end
		end
	end
	# if good_num>0
	#     @info "erro_______________________________________________________________________$good_num _____________________________________"
	# end
	return good_transformation, UpdateNeibors
end
function payload_schedule_to_action_schedule_inf(connectivity, payload_schedule_sim)
	N = size(connectivity, 1)
	state = zeros(Int8, N, N, N)
	for node_id ∈ 1:N
		state[node_id, 1, node_id] = 1
	end
	action_schedule = []
	for payload_simu in payload_schedule_sim
		if payload_simu["pt"] == "g"
			payload = payload_simu["p"]
			tx_id = payload_simu["t"]
			action = coded_payload_to_action(state, tx_id, payload)
			apply_action!(connectivity, state, action)
			purge_disseminated_contents!(state, get_is_decoded(state))
			action_dict = Dict(
				"at" => "g",
				"a" => action,
			)
			push!(action_schedule, action_dict)
		elseif payload_simu["pt"] == "m"
			payload_s = payload_simu["p"]
			sim_action = []
			index = 1
			tx_list = payload_simu["t"]
			for x in payload_s # 同时传输的负载
				payload = x
				tx_id = tx_list[index]
				action = coded_payload_to_action(state, tx_id, payload)
				# @info action
				if action != false
					apply_action!(connectivity, state, action)
					purge_disseminated_contents!(state, get_is_decoded(state))
					push!(sim_action, action)
				end
				index = index = index + 1
			end

			action_dict = Dict(
				"at" => "m",
				"a" => sim_action,
				# "sim_tx_number" => length(sim_action)
			)
			push!(action_schedule, action_dict)
		end
	end
	return action_schedule
end
