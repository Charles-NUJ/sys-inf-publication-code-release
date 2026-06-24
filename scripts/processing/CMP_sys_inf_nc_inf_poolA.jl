#C:\age-optimal-multisource-flooding\data\exp_raw_sys_inf 

#非启发式
#poolcmp  sys nc ncinf  去除8个 去除原来nc中没有的  3*845 = 2535
# 启发式C:\age-optimal-multisource-flooding\data\exp_raw_sys_inf_e 
#poolcmpe sys 和 sysinf 的启发式，来自pool_good  去除原来nc中没有的  2*845=1690

# 131 216 225 265 422 446 601 649
# 7的 873-8=845 

using DrWatson
quickactivate("age-optimal-multisource-flooding")

using DataFrames, Statistics, PyPlot, Graphs, GraphIO

function pl_eps_comp(graph_sizes)

	PyPlot.matplotlib.rcParams["text.usetex"] = false
	PyPlot.matplotlib.rcParams["font.family"] = "serif"
	PyPlot.matplotlib.rcParams["font.serif"] =
		["Times New Roman", "Times", "DejaVu Serif"]
	PyPlot.matplotlib.rc("pdf", fonttype=42)
	PyPlot.matplotlib.rc("ps", fonttype=42)

	fig, ax = PyPlot.subplots(figsize=(4.15, 3.12))
	data_sys_eps0 = []
	data_sys_eps1 = []
	# data_sys_eps2 = []

	data_nc_eps0 = []
	data_nc_eps1 = []
	# data_nc_eps2 = []

	# data_pg_eps0 = []
	# data_pg_eps1 = []
	# data_pg_eps2 = []


# sys 	nc 		pg
# 1  	 n 		 n  # max_inflight
# false true  true # piggyback_enable
	for (i, graph_size) in enumerate(graph_sizes)
		# df = collect_results(datadir("exp_raw_eps/$top_type/$graph_size"))
		df = collect_results("C:/age-optimal-multisource-flooding/data/exp_raw_sys_inf/poolcmpinf_e/$graph_size")
		df = groupby(df, :max_inf, sort = true)

		sys = df[1]
		inf = df[2]

		# @info sys
		# @info inf

		df2 = collect_results("C:/age-optimal-multisource-flooding/data/exp_raw_sys_inf/poolnc3e/$graph_size") # sys nc ncinf
		df2 = groupby(df2, :max_inf, sort = true)
		nc =  df2[1]
		ncinf = df2[2]



		aoi_sys = [direct[1][1] for direct in sys.aoi_from_direct]
		aoi_inf = [direct[1][1] for direct in inf.aoi_from_direct]

		aoi_nc =  [direct[1][1] for direct in nc.aoi_from_direct]
		aoi_ncinf =  [direct[1][1] for direct in ncinf.aoi_from_direct]


		# @info aoi_sys

		for i in aoi_sys
			push!(data_sys_eps0, i)
		end


		for i in aoi_inf
			push!(data_sys_eps1, i)
		end


		for i in aoi_nc
			push!(data_nc_eps0, i)
		end


		for i in aoi_ncinf
			push!(data_nc_eps1, i)
		end

	end

	sort!(data_sys_eps0)
	@info "sys mean  $(mean(data_sys_eps0))"

	sort!(data_sys_eps1)
	@info "sys inf  $(mean(data_sys_eps1))"

	sort!(data_nc_eps0)
	@info "nc   $(mean(data_nc_eps0))"
	sort!(data_nc_eps1)
	@info "nc inf  $(mean(data_nc_eps1))"


	# [ Info: sys mean  10.111469221835073
	# [ Info: sys inf  9.721385017421591
	# [ Info: nc   9.541906213704992
	# [ Info: nc inf  9.220204219899342
	# xs = 1:1:27#984
	# xs = 1:1:6#984
	xs = 1:1:984
	
	# # xs = 4:1:10
	# # data = [data_sys_eps0 data_sys_eps1 data_sys_eps2];
	# # data = hcat(data_sys_eps0, data_sys_eps1, data_sys_eps2,data_pg_eps0, data_pg_eps1, data_pg_eps2)
	# # data = hcat(data_sys_eps0, data_sys_eps1, data_sys_eps2, data_nc_eps0, data_nc_eps1, data_nc_eps2, data_pg_eps0, data_pg_eps1, data_pg_eps2)
	data = hcat(data_sys_eps0, data_sys_eps1,data_nc_eps0, data_nc_eps1)
	# # title_top = ""

	# @info data[:,1]
	# @info data[:,2]
	# @info data[:,3]
	# @info data[:,4]
	Marker_Size = 2
	Line_Width=2
	# # Plot the data, specifying colors individually for each line

	plot(xs, data[:, 1], label = "sys",  linestyle = "-", color = "black", linewidth = Line_Width)
	plot(xs, data[:, 2], label = "sys-inf", linestyle = "-", color = "blue", linewidth = Line_Width)

	# plot(xs, data[:, 1], label = "sys", marker = "o", markersize = Marker_Size, linestyle = "-", color = "black", linewidth = Line_Width)
	# plot(xs, data[:, 2], label = "sys-inf", marker = "o", markersize = Marker_Size, linestyle = "--", color = "blue", linewidth = Line_Width)
	# plot(xs, data[:, 3], label = "\$\\epsilon=0.5\$, sys", marker = "o", markersize = Marker_Size, linestyle = "--", color = "red", linewidth = 1)


	plot(xs, data[:, 3], label = "nc",  linestyle = "-", color = "green", linewidth = Line_Width)
	plot(xs, data[:, 4], label = "nc-inf", linestyle = "-", color = "red", linewidth = Line_Width)




	# plot(xs, data[:, 3], label = "nc", marker = "s", markersize = Marker_Size, linestyle = "-", color = "gray", linewidth = Line_Width)
	# plot(xs, data[:, 4], label = "nc-inf", marker = "s", markersize = Marker_Size, linestyle = "--", color = "red", linewidth = Line_Width)
	# plot(xs, data[:, 6], label = "\$\\epsilon=0.5\$, nc, \$ L_{max}=2\$", marker = "s", markersize = Marker_Size, linestyle = "--", color = "red", linewidth = 1)

	# plot(xs, data[:, 7], label = "\$\\epsilon=0\$, pg, \$L_{max}=2\$", marker = "^", markersize = Marker_Size, linestyle = "-", color = "black", linewidth = 1) #>
	# plot(xs, data[:, 8], label = "\$\\epsilon=0.25\$, pg, \$L_{max}=2\$", marker = "^", markersize = Marker_Size, linestyle = "--", color = "blue", linewidth = 1) #>
	# plot(xs, data[:, 9], label = "\$\\epsilon=0.5\$, pg, \$L_{max}=2\$", marker = "^", markersize = Marker_Size, linestyle = "--", color = "red", linewidth = 1)



	# if top_type == "line"
	# 	title_top = "Line"
	# 	yticks(collect(0:30:270), fontsize = 10)
	# elseif top_type == "circle"
	# 	title_top = "Circle"
	# 	yticks(collect(0:30:270), fontsize = 10)
	# elseif top_type == "full"
	# 	title_top = "Full-Connected"
	# 	   yticks(collect(0:10:70), fontsize = 10)
	# end



	# # Set plot labels and title
	# xlabel("Network Size \$ N\$ in $title_top Network", fontsize = 14)
	xlabel("sorted topology index", fontsize = 12)
	
	# # xlabel("\$ N\$", fontsize = 12,fontfamily  = "Times New Roman")

	# xticks(collect(3:1:15), fontsize = 10, fontfamily = "Times New Roman")



	yticks(fontsize = 12)
	# # ylabel("Average of Age of Information (AoI)")
	# # ylabel("Average Age of Information (\$\\Delta_{avg} \$) ", fontsize = 12)
	# ylabel("Average Age of Information (\$\\Delta_{avg} \$) ", fontsize = 12)
	ylabel("Average Age of Information (\$\\bar{\\Delta}\$) ", fontsize = 12)
	# # ylabel(" \$\\Delta_{avg} \$ ", fontsize = 12,fontfamily  = "Times New Roman")
    # ylabel("Average Age of Information (\$\\Delta_{avg} \$) ", fontsize = 12,fontfamily  = "Times New Roman")
	# # title("Age of Information vs Graph Size for Different Epsilon Values")
	# # title("\$\\Delta_{avg} \$ vs Graph Size for Different \$ \\epsilon in $title_top Topology Network \$")

	# # Show the legend
	legend(loc = "upper left", fontsize = 12)
	# # yticks(collect(-25:25:75), fontsize = 12)
	savefig(plotsdir("aoi-CMP_sys_inf_nc_inf_pool2A.png"), bbox_inches = "tight")
	savefig(plotsdir("aoi-CMP_sys_inf_nc_inf_pool2A.pdf"), bbox_inches = "tight")
	savefig(plotsdir("aoi-CMP_sys_inf_nc_inf_pool2A.eps"), bbox_inches = "tight")
end
# graph_sizes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
# graph_sizes = [4,5,6,7,8,9,10]#,7,8,9,10]
# graph_sizes =[3, 5, 7, 9, 11, 13, 15]

# top_type = "line"
# pl_eps_comp(graph_sizes, top_type)

graph_sizes=[4,5,6,7]
# top_type = "circle"
# pl_eps_comp(graph_sizes, top_type)
# top_type = "full"
pl_eps_comp(graph_sizes)
# # include("scripts/processing/EPS_sys_nc_line_circle.jl")
# # include("scripts/processing/EPS_sys_nc_line_circle-aoi-by-dist.jl")
