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

using Trapz, LinearAlgebra
include("$(srcdir())/utilities/state_utilities_src_inf_v3.jl")

function update_first_decoded!(first_decoding, is_decoded, action_id)
	N = size(is_decoded, 1)
	for node_id ∈ 1:N, content_id ∈ 1:N
		if is_decoded[content_id, node_id] == 1
			first_decoding[content_id, node_id] = min(first_decoding[content_id, node_id], action_id + 1)
		end
	end
	return first_decoding
end

function get_aoi(first_encodings, first_decodings, T)
	N = size(first_decodings, 1)
	aoi = zeros(Int16, 2 * T + 1, N, N)

	# Place decoding event aois
	for content_id ∈ 1:N, receiver_id ∈ 1:N
		if content_id != receiver_id
			decoded_at = first_decodings[content_id, receiver_id]
			encoded_at = first_encodings[content_id]
			aoi[decoded_at, content_id, receiver_id] = decoded_at - encoded_at
			aoi[decoded_at+T, content_id, receiver_id] = decoded_at - encoded_at
		end
	end

	# Increment aoi by one for every timestep from decoding events
	for step_id ∈ 2:2*T, content_id ∈ 1:N, receiver_id ∈ 1:N
		# If zero, then the value was never set
		if aoi[step_id, content_id, receiver_id] == 0
			aoi[step_id, content_id, receiver_id] = aoi[step_id-1, content_id, receiver_id] + 1
		end
	end

	# Interspace aoi trace with two values per timestep
	double_valued_aoi = zeros(2 * T, N, N)
	double_valued_time = zeros(Float64, 2 * T)
	for step ∈ 1:T, node_id ∈ 1:N, content_id ∈ 1:N
		double_valued_aoi[(step-1)*2+1, content_id, node_id] = aoi[T+step, content_id, node_id]
		double_valued_aoi[step*2, content_id, node_id]       = aoi[T+step, content_id, node_id] + 1
		double_valued_time[(step-1)*2+1]                     = step
		double_valued_time[step*2]                           = step + 1
	end

	return double_valued_aoi, double_valued_time
end

function calculate_avg_aoi(first_encodings, first_decodings, T)
	double_valued_aoi, double_valued_time = get_aoi(first_encodings, first_decodings, T)
	N = length(first_encodings)
	aoi_summer = 0.0
	for node_id ∈ 1:N, content_id ∈ 1:N
		if node_id != content_id
			aoi_summer += trapz((double_valued_time), double_valued_aoi[:, content_id, node_id]) / T
		end
	end
	return aoi_summer / (N * (N - 1))
end

function get_coding_events(connectivity, schedule)
	N = size(connectivity, 1)
	T = size(schedule, 1)
	# content_id x node_id
	first_decodings = ones(Int16, N, N) .* 32767
	first_decodings[I(N)] .= 0
	first_encodings = ones(Int16, N) .* 32767
	state = zeros(Int16, N, N, N)
	for node_id ∈ 1:N
		state[node_id, 1, node_id] = 1
	end

	for (step, action) in enumerate(schedule)
		apply_action!(connectivity, state, action)
		tx_id = action[1]
		payload = action_to_coded_payload(state, action)
		if payload[tx_id] == 1
			first_encodings[tx_id] = min(first_encodings[tx_id], step)
		end
		is_decoded = get_is_decoded(state)
		update_first_decoded!(first_decodings, is_decoded, step)
		purge_disseminated_contents!(state, is_decoded)
	end

	return first_encodings, first_decodings
end


function get_aoi_components(connectivity, schedule)
	N = size(connectivity, 1)
	T = size(schedule, 1)
	first_encodings, first_decodings = get_coding_events(connectivity, schedule)
	avg_aoi = calculate_avg_aoi(first_encodings, first_decodings, T)
	dissemination_delays = zeros(Int16, N, N)
	for content_id ∈ 1:N
		dissemination_delays[content_id, :] = first_decodings[content_id, :] .- first_encodings[content_id]
	end
	dissemination_delays[I(N)] .= 0

	aoi_components = (avg_aoi, T / 2, sum(dissemination_delays) / N / (N - 1))
	return aoi_components
end


# 应用动作并更新状态，记录首次编码和解码时间
function get_aoi_components_direct(connectivity, schedule, first_encodings, first_decodings)
	N = size(connectivity, 1)
	T = size(schedule, 1)
	# first_encodings, first_decodings = get_coding_events(connectivity, schedule) # 计算第一次发送 和 第一次相互接收到的时间
	avg_aoi = calculate_avg_aoi(first_encodings, first_decodings, T)
	dissemination_delays = zeros(Int16, N, N)
	for content_id ∈ 1:N
		dissemination_delays[content_id, :] = first_decodings[content_id, :] .- first_encodings[content_id]
	end
	dissemination_delays[I(N)] .= 0
	aoi_components = (avg_aoi, T / 2, sum(dissemination_delays) / N / (N - 1))
	return aoi_components
end


# 表示正在传，但是还没有被所有节点收到的信息
function get_current_inflight_data(step, first_encodings, first_decodings)
	N = length(first_encodings)
	inflight = []

	# all_decoded_judge = []


	for c_id ∈ 1:N
		was_txed = (step > first_encodings[c_id]) #一个节点内容发出
		all_decoded = true


		for n_id ∈ 1:N
			if c_id != n_id
				if step <= first_decodings[c_id, n_id]# 默认信息是127 因此，在任何一个其余节点上没有被收到，都会被标记
					all_decoded = false
					# push!(all_decoded_judge,1)
				end
			end
		end

		# if length(all_decoded_judge) > 0
		# if sum(all_decoded_judge)  >= N-1

		if was_txed && !all_decoded # 信息已经发出，而且在有的节点上没有被收到
			append!(inflight, c_id)
		end
		# end
		# end


	end
	return inflight
end


# 表示之前各个时隙中最大的信息传输量（用于 /process/inflight.jl 最终结果统计中 inflight 最大值比例的寻找）
function get_max_inflight_data(first_encodings, first_decodings)
	T = maximum(first_decodings)# 二维矩阵中最大的数,代表最大的接收到的时间
	N = length(first_encodings)
	max_inflight = 0
	for t ∈ 1:T # 之前的所有传输时隙中
		num_inflight = 0
		for c_id ∈ 1:N
			was_txed = (t > first_encodings[c_id])
			all_decoded = true
			for n_id ∈ 1:N
				if c_id != n_id
					if t <= first_decodings[c_id, n_id]
						all_decoded = false
					end
				end
			end
			if was_txed && !all_decoded
				num_inflight += 1
			end
		end
		# 各个时隙中最大的同时传输量
		max_inflight = max(max_inflight, num_inflight)
	end
	return max_inflight
end
