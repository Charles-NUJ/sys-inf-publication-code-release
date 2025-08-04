using Graphs, GraphIO, Distributed, DataFrames, Statistics
# include("$(srcdir())/solvers/a_star.jl")
# include("$(srcdir())/solvers/a_star_src.jl")  # 1
# include("$(srcdir())/solvers/a_star_src_inf.jl") # 2
# include("$(srcdir())/solvers/incremental_a_star.jl") # 3 
# include("$(srcdir())/solvers/incremental_a_star_inf.jl") # 4
# include("$(srcdir())/solvers/incremental_a_star_inf_v2.jl") # 4
include("$(srcdir())/solvers/incremental_a_star_inf_v3.jl") # 4
# include("$(srcdir())/solvers/incremental_a_star_inf_v4.jl") # 4
include("$(srcdir())/solvers/nc_capacities.jl")
include("$(srcdir())/utilities/topology_utilities.jl")

# avg_outdegree = mean(outdegree(graph)


function get_aoi_by_trapezoid(schedule, first_encodings, first_decodings)
	T = length(schedule)
	N = length(first_encodings)
	dissemination_delays = zeros(Int16, N, N)

	for content_id ∈ 1:N
		dissemination_delays[content_id, :] = first_decodings[content_id, :] .- first_encodings[content_id]
		# solution_knowledge["dissemination_delays"] = first_decodings .- first_encodings 在引入不编码时候的界的时候是这么计算的
	end
	dissemination_delays[I(N)] .= 0
	@info dissemination_delays

	S_trapezoid = []
	for i ∈ 1:N, j ∈ 1:N
		if i != j
			D = dissemination_delays[i, j]
			push!(S_trapezoid, (T + 2D) / 2)
		end
	end

	S_avg = S_trapezoid # 在一个周期内的平均面积

	aoi = sum(S_avg) / N / (N - 1)
	# @info aoi
	return aoi
end

function get_aoi_components_direct(connectivity, schedule, first_encodings, first_decodings)
	N = size(connectivity, 1)
	T = size(schedule, 1)
	# first_encodings, first_decodings = get_coding_events(connectivity, schedule) # 计算第一次发送 和 第一次相互接收到的时间
	avg_aoi = calculate_avg_aoi_direct(first_encodings, first_decodings, T)
	dissemination_delays = zeros(Int16, N, N)
	for content_id ∈ 1:N
		dissemination_delays[content_id, :] = first_decodings[content_id, :] .- first_encodings[content_id]
	end
	dissemination_delays[I(N)] .= 0
	aoi_components = (avg_aoi, T / 2, sum(dissemination_delays) / N / (N - 1))
	return aoi_components
end
function get_aoi_direct(first_encodings, first_decodings, T)
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
function calculate_avg_aoi_direct(first_encodings, first_decodings, T)
	double_valued_aoi, double_valued_time = get_aoi_direct(first_encodings, first_decodings, T)
	N = length(first_encodings)
	aoi_summer = 0.0
	for node_id ∈ 1:N, content_id ∈ 1:N
		if node_id != content_id
			aoi_summer += trapz((double_valued_time), double_valued_aoi[:, content_id, node_id]) / T
		end
	end
	return aoi_summer / (N * (N - 1))
end
function new_get_aoi(graph_size, schedule, first_encodings, first_decodings)
	N = graph_size
	T = size(schedule, 1)
	dissemination_delays = zeros(Int16, N, N)
	for content_id ∈ 1:N
		dissemination_delays[content_id, :] = first_decodings[content_id, :] .- first_encodings[content_id]
	end

	dissemination_delays[I(N)] .= 0
	# @info "dissemination_delays $dissemination_delays" 
	return T / 2 + sum(dissemination_delays) / N / (N - 1)
end
function run_experiment(graph_size, config)

	connectivity = []
	graph_type = config["graph_type"]
	Diameter = 1
	# Load experiment config
	if graph_type == "pool"
		graphs = readlines("$(datadir())/exp_res/graph$(graph_size)c$(config["cutoff"]).g6")
		graph = GraphIO.Graph6._g6StringToGraph(graphs[config["graph_id"]])
		connectivity = collect(adjacency_matrix(graph))
	elseif graph_type == "line" || graph_type == "circle"

		N = graph_size
		connectivity = zeros(Int8, N, N)
		for n in 1:N-1
			connectivity[n, n+1] = 1
		end
		for n in 2:N
			connectivity[n, n-1] = 1
		end

		if graph_type == "circle"
			connectivity[1, N] = 1
			connectivity[N, 1] = 1
		end

	elseif graph_type == "euclidean_graph" || graph_type == "erdos_renyi"

		if graph_type == "euclidean_graph"
			graphs = readlines("$(datadir())/exp_res/graph$(graph_size)c$(config["cutoff"]).g6")
		else
			graphs = readlines("$(datadir())/exp_res/renyi_graph$(graph_size)c$(config["cutoff"]).g6")
		end
		graph = GraphIO.Graph6._g6StringToGraph(graphs[config["graph_id"]])
		connectivity = collect(adjacency_matrix(graph))
		# @info connectivity
		# end
	end

	solution_knowledge = Dict([("dissemination_delays" => zeros(graph_size, graph_size)), ("t_star", 0)]) # Lower Bound Initialisation
	# @info "6 文件 graph_analysis 拓扑" connectivity
	# @info "7" solution_knowledge
	# Prepare solver input, with or without network coding
	if config["max_inference_simu"] == 1 # 非干扰的为1

		println("Running 1：Farazi")
		mcds_cardinality, non_leaves, true_leaves, pseudo_leaves, mcds = get_node_types(connectivity)
		solution_knowledge["t_star"] = graph_size * mcds_cardinality + (graph_size - length(non_leaves))


	else # 这里提供了下界
		println("Running 2：Farazi")
		mcds_cardinality, non_leaves, true_leaves, pseudo_leaves, mcds = get_node_types(connectivity)
		solution_knowledge["t_star"] = graph_size * mcds_cardinality + (graph_size - length(non_leaves))


		# 使用的界限 为Widmer论文的界限
		# println("Running 2: Widmer")
		# t_star taken from Widmer paper
		# solution_knowledge["t_star"] = sum(min_schedule_length_nc(connectivity))

		# 这里的schedule 长度不是这里决定的
		# We might need to have information from the non-network coding solution, so try to load it

		using_sys_band = false
		# using_sys_band = true
		if using_sys_band
			sys_config = deepcopy(config)
			sys_config["max_inflight"] = 1
			sys_config["max_inference_simu"] = 1
			sys_config["inf"] = false
			sys_results = missing


			# println("$sys_config")
			# println("$config")
			# !config["heuristic"]

			# if (graph_type == "pool" || graph_type == "line" || graph_type == "circle") || !config["heuristic"]


			# 	if config["heuristic"]

			# 		sys_config["heuristic"] = false

			# 	end

			# 	sys_results = try
			# 		wload(datadir("exp_raw_sys_inf/$(graph_type)/$(graph_size)/$(savename(sys_config,"jld2"))"))
			# 	catch e
			# 		@info "missing ___________________sys dissemination_delays________________________"
			# 		sys_results = missing
			# 	end

			# else
			if graph_type == "euclidean_graph" || graph_type == "erdos_renyi"

				sys_results = try
					wload(datadir("exp_raw_sys_inf_e/$(graph_type)/$(graph_size)/$(savename(sys_config,"jld2"))"))
				catch e
					sys_results = missing
				end

			end

			if !ismissing(sys_results) # 没有就是miss, 有就是 !ismissing; 有不编码的结果话, 获取下界，用来加速收敛
				# first_encodings, first_decodings = get_coding_events(connectivity, sys_results["schedule"])
				@info "using ___________________sys as sys-inf dissemination_delays band ________________________"
				first_encodings = sys_results["first_encodings"]
				first_decodings = sys_results["decoding_levels"]

				# decoding_levels[I(graph_size)] .= 0
				# first_decodings = decoding_levels
				solution_knowledge["dissemination_delays"] = first_decodings .- first_encodings # Tighten lower bounds to improve convergence
				# @info first_decodings
				# @info first_encodings
				solution_knowledge["dissemination_delays"][I(graph_size)] .= 0
				# @info solution_knowledge["dissemination_delays"]
			else
				@info "error——————————————————————————————————————————————————————————————"
				# solution_knowledge["dissemination_delays"]=
			end
		end
	end

	if graph_type == "pool" || graph_type == "euclidean_graph" || graph_type == "erdos_renyi"
		Diameter = diameter(graph)
	else
		graph = SimpleDiGraph(connectivity)
		Diameter = diameter(graph)
	end


	errorFlage =  false

	# print("solution_knowledge ", solution_knowledge, "\n")
	if (config["max_inference_simu"] == 1) && !config["heuristic"] && !config["piggyback_enable"] && !config["inf"]

		println("1 A star")
		# processing_time = @elapsed best_solution = a_star_search(connectivity, config["max_coding_degree"], config["max_coding_degree"], deepcopy(solution_knowledge))

	elseif (config["max_inference_simu"] >= 2) && !config["heuristic"] && !config["piggyback_enable"] && config["inf"]

		println("2 A star inf")
		# processing_time = @elapsed best_solution = a_star_inf_search(connectivity, config["max_coding_degree"], config["max_inflight"], config["max_inference_simu"], deepcopy(solution_knowledge))

	elseif (config["max_inference_simu"] == 1) && config["heuristic"] && !config["piggyback_enable"] && !config["inf"]

		println("3 A star e")
		# processing_time = @elapsed best_solution = incremental_a_star_search(connectivity, config["max_coding_degree"], config["max_coding_degree"], deepcopy(solution_knowledge))

	elseif (config["max_inference_simu"] >= 2) && config["heuristic"] && !config["piggyback_enable"] && config["inf"]

		println("4 A star inf-e")
		set_max_inference_simu = 2


		# Color  =  greedy_color(graph; sort_degree=true, reps = 1) 
		# x = Color.colors
		# colornum =  Color.num_colors
		# @info "colornum" colornum
		# s=counter(x)
		# list_num = []
		# for i in s
		# 	push!(list_num,i[2])
		# end

		# X = maximum(list_num)

		if Diameter <= 2
			set_max_inference_simu = 1
		elseif Diameter <= 5
			set_max_inference_simu = 2
		elseif Diameter <= 8
			set_max_inference_simu = 3
		elseif Diameter <= 11
			set_max_inference_simu = 4
		elseif Diameter <= 14
			set_max_inference_simu = 5
		elseif Diameter <= 17
			set_max_inference_simu = 6
		elseif Diameter <= 20
			set_max_inference_simu = 7
		end
		# set_max_inference_simu = X
		# println("Graph Diameter $Diameter  X = $set_max_inference_simu")
		println("Graph Diameter $Diameter  set_max_inference_simu = $set_max_inference_simu")

		if graph_type == "euclidean_graph" || graph_type == "erdos_renyi" ||graph_type=="line"
			if set_max_inference_simu >= 3
				set_max_inference_simu = 2
				println("====================> set_max_inference_simu $set_max_inference_simu")
				# @info "______________________________________________________________________________________________"
			end
		end

		@info "$(config["max_coding_degree"]), $set_max_inference_simu, $set_max_inference_simu"
		# processing_time = @elapsed best_solution = incremental_a_star_inf_search(connectivity, config["max_coding_degree"], config["max_inflight"], config["max_inference_simu"], deepcopy(solution_knowledge))

		try
			processing_time = @elapsed best_solution = incremental_a_star_inf_search(connectivity, config["max_coding_degree"], set_max_inference_simu, set_max_inference_simu, deepcopy(solution_knowledge))
			errorFlage =  false

			log_result = false
			log_result = true
			if log_result
				# @info "结束 $best_solution"
				schedule = best_solution.schedule
				aoi_from_dist = best_solution.g_score ./ (graph_size * (graph_size - 1))
		
				aoi_components = 0
				payload_schedule = best_solution.payload_schedule
				first_encodings = best_solution.first_encodings
				decoding_levels = best_solution.decoding_levels
				decoding_levels_temp = decoding_levels
				decoding_levels_temp[I(graph_size)] .= 1
				first_decodings = decoding_levels_temp
		
				println("-------- result -------------------")
				# print("first_encodings ", first_encodings, "\n")
				# print("decoding_levels ", decoding_levels, "\n")
				# print("first_decodings ", first_decodings, "\n")
				# print("best_solution ", best_solution, "\n")
				aoi_from_direct = get_aoi_components_direct(connectivity, schedule, first_encodings, first_decodings)
				aoi_by_trapezoid = get_aoi_by_trapezoid(schedule, first_encodings, first_decodings)
				aoi_get_new = new_get_aoi(graph_size, schedule, first_encodings, first_decodings)
				# aoi_components = get_aoi_components(connectivity, schedule)
				max_piggy = 0.0
				# print(aoi_components, "\n")
				# Collecting Result Dataframe
				graph_id = config["graph_id"]
				# avg_outdegree = mean(outdegree(graph))
				# graph = graphs[config["graph_id"]]
				cutoff = config["cutoff"]
				t = length(schedule)
				t_star = solution_knowledge["t_star"]
				max_coding_degree = config["max_coding_degree"]
				max_inflight = config["max_inflight"]
				heuristic = config["heuristic"]
		
				piggyback_enable = config["piggyback_enable"]
				inf_enable = config["inf"]
				eps = config["eps"]
				max_inf = config["max_inference_simu"]
				max_inf_total = 1
				if config["inf"]
					for schedule_simu in schedule
						if schedule_simu["at"] == "m"
							# tx_sim_number = schedule_simu["sim_tx_number"]
							tx_sim_number = length(schedule_simu["a"])
							max_inf_total = max(max_inf_total, tx_sim_number)
						end
					end
				end
				result =
					@strdict graph_id schedule t t_star aoi_by_trapezoid aoi_from_dist aoi_from_direct aoi_components max_coding_degree heuristic max_inflight max_inf max_piggy processing_time cutoff payload_schedule decoding_levels first_encodings eps piggyback_enable inf_enable max_inf_total Diameter aoi_get_new
		
				if heuristic
					wsave(datadir("exp_raw_sys_inf_e/$(config["graph_type"])/$(graph_size)", savename(config, "jld2")), result)
				else
					wsave(datadir("exp_raw_sys_inf/$(config["graph_type"])/$(graph_size)", savename(config, "jld2")), result)
				end
		
				# @info result
				# @info schedule
				# @info payload_schedule
				@info "Results \n ##EPS: $eps \n ## Graph Size: $graph_size \n ## Graph ID: $graph_id  \n ## Piggyback: $(piggyback_enable)  \n ## max_coding_degree: $(max_coding_degree) \n ## max_inflight: $(max_inflight) \n ## max_inf: $(max_inf)	 \n ## max_piggy: $(max_piggy)  \n ## heuristic: $(heuristic) \n ## inf: $(inf_enable) \n ## max_inf_total: $(max_inf_total) \n ## Schedule Length: $(t) \n ## Avg. AoI: $(aoi_components[1]) \n ## aoi_by_trapezoid $aoi_by_trapezoid \n ## AoI aoi_from_dist: $(aoi_from_dist) \n ## AoI aoi_from_direct: $(aoi_from_direct) \n ## AoI aoi_get_new: $(aoi_get_new) \n ## Processing Time $processing_time \n ## schedule[0] $(schedule[1])"
				println("")
			end

		catch e
			errorFlage =  true
		end
		# else
		# 	 # 把sys跑的复制过来即可
		# 	sys_config = deepcopy(config)
		# 	sys_config["max_inflight"] = 1
		# 	sys_config["max_inference_simu"] = 1
		# 	sys_config["inf"] = false

		# 	dir_src = datadir("exp_raw_sys_inf_e/$(config["graph_type"])/$(graph_size)", savename(sys_config, "jld2"))
		# 	dir_dst = datadir("exp_raw_sys_inf_e/$(config["graph_type"])/$(graph_size)", savename(config, "jld2"))
		# 	# dir_src = datadir("", savename(sys_config, "jld2"))
		# 	# dir_dst = datadir("", savename(config, "jld2"))
		# 	cp(dir_src,dir_dst,force=true)
		# end

	else
		@info "error config"
	end

	if !errorFlage
	
	end
	return nothing
end
