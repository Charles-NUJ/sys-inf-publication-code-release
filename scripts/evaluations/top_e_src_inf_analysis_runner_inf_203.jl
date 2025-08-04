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
# include("scripts/evaluations/top_e_src_inf_analysis_runner_inf.jl")
# include("scripts/evaluations/Inf_analysis_runner.jl")
using DrWatson
quickactivate("age-optimal-multisource-flooding")

using Graphs, GraphIO, Distributed
include("$(srcdir())/utilities/topology_utilities.jl")
include("$(scriptsdir())/evaluations/top_e_src_inf_graph_analysis_inf.jl")
using Dates: Dates #println("$(Dates.now(Dates.UTC))")
# 来自数据库的并发传输, 根据 graph_type 进行区分不同的拓扑, pool, euclidean_graph,erdos_renyi,line,circle
# 应对sys,sys-e, sys-inf,sys-inf-e
# pool 对应 sys 和 sys-inf 节点数 3-7 -- 绘制增益矩形图, 以及 sys-inf的 max-inflaght 和 max_inference_simu 的表格, 这里存在一个问题是 max大于1时候,计算编码代价还在代码里, 需要对 sys-inf a-satr-inf修改
# 以及重新跑一些图

# 两种随机拓扑对应于 sys-e 和 sys-inf-e 节点数为 8:2:20-- 绘制热度增益图

# line 和 circle 图用于初步的结果呈现
# 3-7 用 sys 和 sys-inf
# 8:2:20 用sys-e 和 sys-inf-e  这里同样要注意, 对于  max大于1时候, 编码代价还在代码里需要重新跑图


# 在Jula-Run 下面有绘图的代码

function run(graph_sizes, force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable, num_random_graphs = 1)
	for graph_type in graph_types
		# Setup
		for graph_size in graph_sizes # Number of nodes
			piggyback_enable = [false]
			max_coding_degrees = [1] # 最大编码度, 表示节点本地一次可以传输的包中, 含有的信息量
			max_inflights = [1] # 最大飞行量, 表示网络中正在传输但是还没有所有节点都收到的信息个数
			max_inference_simu = [1] # 默认支持最多并发，表示网络中允许的同时传输的节点数

			# Runner Config
			# force_calculation = false # Overwrite any existing results

			if inf_enable[1] == true
				if graph_type == "line"
					max_inflights = [2]
					max_inference_simu = [2]
					# if graph_size <= 3
					# 	max_inflights = [2, 3]
					# 	max_inference_simu = [2, 3]
					# elseif graph_size <= 6
					# 	max_inflights = [2]
					# 	max_inference_simu = [2]
					# elseif graph_size <= 9
					# 	max_inflights = [2]
					# 	max_inference_simu = [2]
					# elseif graph_size <= 12
					# 	max_inflights = [2]#,34] 
					# 	max_inference_simu = [2]#,3]#4]
					# elseif graph_size <= 15
					# 	max_inflights = [2]#,3]#5] 
					# 	max_inference_simu = [2]#,3]#5]
					# elseif graph_size <= 18
					# 	max_inflights = [2]#,3]#6] 
					# 	max_inference_simu = [2]#,3]#6]
					# elseif graph_size <= 21
					# 	max_inflights = [2]#,3]#7] 
					# 	max_inference_simu = [2]#,3]#7]
					# end
				elseif graph_type == "circle"
					if graph_size <= 5
						max_inflights = [2, 3]
						max_inference_simu = [2, 3]
					elseif graph_size <= 8
						max_inflights = [2, 3]
						max_inference_simu = [2, 3]
					elseif graph_size <= 11
						max_inflights = [2, 3]
						max_inference_simu = [2, 3]
					elseif graph_size <= 14
						max_inflights = [2]#,34] 
						max_inference_simu = [2]#,3]#4]
					elseif graph_size <= 17
						max_inflights = [2]#,3,4]#5] 
						max_inference_simu = [2]#,3,4]#5]
					elseif graph_size <= 20
						max_inflights = [2]#,3,4,5]#6] 
						max_inference_simu = [2]#,3,4,5]#6]
					elseif graph_size <= 23
						max_inflights = [2]#,3,4,5]#7] 
						max_inference_simu = [2]#,3,4,5]#7]
					end
				elseif graph_type == "pool" # 4  5 6 7
					max_inflights = [graph_size]
					max_inference_simu = [graph_size]
					# if graph_size <=3
					# 	max_inflights = [1] 
					# 	max_inference_simu = [1] 
					# else
					if graph_size <= 6
						max_inference_simu = [2]
					elseif graph_size <= 7
						# max_inflights = [2,3] 
						# max_inference_simu = [2,3]
						max_inference_simu = [3]
					end
				elseif graph_type == "euclidean_graph" || graph_type == "erdos_renyi"
					max_inflights = [2]
					max_inference_simu = [2]
				end
			elseif inf_enable[1] == false
				max_inflights = [1]
				max_inference_simu = [1]
			end

			# heuristic = [true] # Enabled: divide-and-conquer(分而治之) method using the max-inflight heuristic, use heuristic for graph_size > 7
			# generate_random_topologies = false # Exhaustive connected graph sets available fir 3 < until graph_size < 8, use random topologies for larger networks 

			# Random Topologies
			# graph_type = "euclidean_graph" # ["euclidean_graph", "erdos_renyi"]
			# num_random_graphs = num_random_graph # 30# Number of random graphs to produced

			## -- ##
			cutoffs = [""]
			num_graphs = num_random_graphs

			if graph_type == "euclidean_graph"
				cutoffs = collect(0.2:0.05:1.0)
			elseif graph_type == "erdos_renyi"
				###########################设置为0.2开始，设置为0.25是为了分到MSI上#####################################################################################
				cutoffs = collect(0.2:0.05:0.95) # Graph property: cut-off distance (Eucliden Graphs) or density (Erdos Renyi Graphs)
				# cutoffs = collect(0.3:0.05:0.95) # Graph property: cut-off distance (Eucliden Graphs) or density (Erdos Renyi Graphs)
				if graph_size == 8
					cutoffs = collect(0.25:0.05:0.95) # Graph property: cut-off distance (Eucliden Graphs) or density (Erdos Renyi Graphs)
				end
			end


			if graph_type == "euclidean_graph" || graph_type == "erdos_renyi"

				if generate_random_topologies[1]  #生成随机拓扑

					data_source = "pool"

					# num_graphs = num_random_graphs
					for cutoff in cutoffs
						seed = 1
						graphs = []
						for _ in 1:num_random_graphs
							while true # Generate a new graph until it is connected
								seed += 1
								if graph_type == "euclidean_graph"
									graph = euclidean_graph(graph_size, 2; cutoff = cutoff, seed = seed)[1]
								elseif graph_type == "erdos_renyi"
									graph = erdos_renyi(graph_size, trunc(Int, cutoff * (graph_size * (graph_size - 1)) / 2), seed = seed)
								else
									println("Wrong graph type selected, available types are: [euclidean_graph, erdos_renyi]")
								end

								if is_connected(graph)
									push!(graphs, GraphIO.Graph6._graphToG6String(graph))
									break
								end
							end
						end

						# for g in graphs
						#     println("$g => $(g[11:end])")
						# end
						# Save graphs to data set
						open("$(datadir())/exp_res/$(graph_type)_$(graph_size)c$(cutoff).g6", "w") do f
							[println(f, g[11:end]) for g in graphs]
						end
					end

					# elseif !generate_random_topologies[1] #载入随机拓扑
					# 	num_graphs = num_random_graphs
				end

			elseif graph_type == "pool"
				num_graphs = size(readlines("$(datadir())/exp_res/graph$(graph_size)c.g6"), 1) # Only checks how many graphs 

				# elseif  graph_type == "line" ||  graph_type == "circle"
				#     num_graphs = 1
			end


			# Experiment Configs
			graph_ids = collect(1:num_graphs)
			# graph_ids = collect(20:num_graphs)
			# @info graph_ids

			# Run (configs are saved so that computation can be manually restarted or passed to a batch system)


			eps = [0.0]
			graph_type_temp = [graph_type]
			graph_size_temp = [graph_size]

			params = dict_list(
				Dict(
					:"graph_type" => graph_type_temp,
					:"graph_size" => graph_size_temp,
					:"graph_id" => graph_ids,
					:"heuristic" => heuristic,
					:"cutoff" => cutoffs,
					:"max_coding_degree" => max_coding_degrees,
					:"max_inflight" => max_inflights,
					:"max_inference_simu" => max_inference_simu,
					:"piggyback_enable" => piggyback_enable,
					:"inf" => inf_enable,
					:"eps" => eps,
					:"regenerate_random" => generate_random_topologies,
				),
			)


			enable_run = true
			if enable_run
				if heuristic[1] == false
					tmp_dir = "tmp_inf/$(graph_type)"
					res = tmpsave(params, projectdir(tmp_dir))
					for r in res
						config = load(projectdir(tmp_dir, r), "params")
						if !isfile(datadir("exp_raw_sys_inf/$(graph_type)/$(graph_size)", savename(config, "jld2"))) || force_calculation
							println("Running exp_raw_sys_inf $config")
							run_experiment(graph_size, config)
						else
							println("Having done $config")
						end

					end
				elseif heuristic[1] == true
					tmp_dir = "tmp_inf_e/$(graph_type)"
					res = tmpsave(params, projectdir(tmp_dir))
					# len = length(res) 1-480
					for r in res#[1:50]
						config = load(projectdir(tmp_dir, r), "params")
						if config["cutoff"]<=0.7 # && (config["graph_id"]==29)# && config["graph_id"]<=5
							if !isfile(datadir("exp_raw_sys_inf_e/$(graph_type)/$(graph_size)", savename(config, "jld2"))) || force_calculation
								println("Running exp_raw_sys_inf_e $config")
								println("$(Dates.now())")
								run_experiment(graph_size, config) #@time 
							else
								println("Having done $config")
							end
						end

					end
				end

			end
		end
	end
end

# include("scripts/evaluations/Inf_analysis_runner.jl") 
# include("scripts/evaluations/top_e_src_inf_analysis_runner_inf.jl")


# graph_sizes = [3,4,5,6,7,8,9,10,11,12,13,14,15,16]#,17,18,19,20,21] #
# heuristic = [false]
generate_random_topologies = [false]
force_calculation = true
# graph_sizes = [9]
# graph_types = ["circle"]
# graph_types = ["line"]
# inf_enable = [false] 
# run(graph_sizes,force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable) # 6.6



# graph_types = ["line"]
# graph_sizes = [11,12,13,14,15,16]
# # # # # # # # # # inf_enable = [false] 
# graph_sizes = [4]#,4,5]
# inf_enable = [true]
# heuristic = [true]
# force_calculation = true
# run(graph_sizes,force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable) # 6.5


# # graph_types = ["circle"]
# # graph_sizes = [3,4,5,6,7,8,9,10,11]
# # graph_sizes = [9]
# # # # inf_enable = [false]
# # inf_enable = [true]
# # run(graph_sizes,force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable)


# ###################### pool ###################
# heuristic = [true]

# graph_sizes = [3, 4, 5, 6, 7]
# # graph_sizes = [7]

# graph_types = ["pool"]
# force_calculation = false
# inf_enable = [true]
# run(graph_sizes, force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable)


############################################### e ####################################################################################
graph_size = [8, 10, 12, 14, 16, 18, 20]
heuristic = [true]

##
generate_random_topologies = [false]
num_random_graphs = 30
heuristic = [true]
force_calculation = true
inf_enable = [true]

graph_sizes = [10]#,16]#,18,20]  v1跑了 8 10 12 和 14的一点点 先当作结果
# graph_types = ["euclidean_graph"]
# run(graph_sizes, force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable, num_random_graphs)

graph_types = ["erdos_renyi"]
run(graph_sizes, force_calculation, graph_types, generate_random_topologies, heuristic, inf_enable, num_random_graphs)
