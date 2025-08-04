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

function get_is_decoded(state)
	N = size(state, 3)
	is_decoded = zeros(Int8, N, N)
	for node_id ∈ 1:N
		for row_id ∈ 1:size(state, 2)
			row_sum = sum(state[:, row_id, node_id])
			if row_sum == 1
			    content_id = findfirst(state[:,row_id,node_id] .==1)
			    is_decoded[content_id, node_id] = 1
			# if row_sum > 0
			# 	is_recived = findall(state[:, row_id, node_id] .> 0)   # 找到有内容的号
			# 	for index in is_recived
			# 		if index > 0
			# 			is_decoded[index, node_id] = 1
			# 		end
			# 	end
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
				# is_recived = findall(state[:, row_id, node_id] .> 0)
				# for index in is_recived
				# 	if index > 0
				# 		is_decoded[content_map[index], node_id] = 1
				# 	end
				# end
			if row_sum == 1
				    content_id = findfirst(state[:,row_id,node_id] .==1)
				    is_decoded[content_map[content_id], node_id] = 1

			elseif row_sum == 0
				break
			end
		end
	end
	return is_decoded
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
				return (tx_id, action)
			end
		end
	end
	return missing
end


function action_to_coded_payload(state, action)
	N = size(state, 3)
	payload = zeros(Int8, size(state, 1))

	tx_id = action[1]
	coding_action = action[2]

	for coding_row ∈ 1:N
		if is_row_encoded(coding_action, coding_row)
			payload .+= state[:, coding_row, tx_id]
		end
	end

	return payload
end

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
	@inbounds for rx_id in connections
		zero_row = findfirst(sum(state[:, :, rx_id], dims = 1) .== 0)
		if !isnothing(zero_row)
			# println("in judge")
			zero_row = zero_row[2]
			for c_d ∈ 1:N
				@inbounds state[c_d, zero_row, rx_id] = (state[c_d, zero_row, rx_id] + payload[c_d]) %2
			end
			bin_mat_rref_line!(state, rx_id, zero_row)
			if sum(state[:, zero_row, rx_id]) > 0
				good_transformation = true
				# good_num = good_num+1
				# println("in judge $good_num")
			end
		end
	end

	# if good_num>0
	#     @info "erro_______________________________________________________________________$good_num _____________________________________"
	# end

	return good_transformation, UpdateNeibors
end

function valid_actions(all_actions, state, max_coding_degree, inactive_nodes = [])
	node_ranks = get_node_ranks(state.decoders)
	step_actions = [action for action in all_actions if (ndigits(action[2], base = 2) <= node_ranks[action[1]]) && (count_ones(action[2]) <= max_coding_degree)]
	# Remove actions of inactive node
	for inactive_node in inactive_nodes
		step_actions = [action for action in step_actions if action[1] != inactive_node]
	end
	return step_actions
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
					if state[content_id, row_id, node_id] == 1
						state[:, row_id, node_id] .= 0
						bin_mat_rref!(state, node_id)
						break
					end
				end
			end
		# end
	end
	return nothing
end

function payload_schedule_to_action_schedule(connectivity, payload_schedule)
	N = size(connectivity, 1)
	state = zeros(Int8, N, N, N)
	for node_id ∈ 1:N
		state[node_id, 1, node_id] = 1
	end
	action_schedule = []
	for payload in payload_schedule
		action = coded_payload_to_action(state, payload[1], payload[2])
		apply_action!(connectivity, state, action)
		purge_disseminated_contents!(state, get_is_decoded(state))
		push!(action_schedule, action)
	end
	return action_schedule
end



