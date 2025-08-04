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

using DrWatson
quickactivate("age-optimal-multisource-flooding")

using DataFrames, Statistics, PyPlot, Graphs, GraphIO, Colors

function rounds(value, step)
    return round(round(value / step) * step,digits=2);
end

# data\exp_raw_e_new\euclidean_graph

function pl_graphsize_degree_gain(graph_sizes, group_by, subresult)
    x_values = []
    y_values = []
    z_values = []
    sample_values =  []

    ## Set up LaTeX fonts
    plt.rcParams["text.usetex"] = true
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serifx"] =  ["Computer Modern Roman"]
    plt.rcParams["font.size"] = 30

    if subresult == "geometric_graphs"
        subresult_path = "euc" # 放了sys 和 nc ecu
        title_string = "Random Geometric Graphs"
        save_title_string ="geometric_graphs"
    else
        subresult_path = "renyi" # 放了sys 和 nc renyi
        title_string = "Erdős-Rényi Graphs"
        save_title_string ="erod_renyi_graphs"
    end

    if group_by == "density"
        bins = 1.0:-0.05:0.2
        # bins = 0.7:-0.05:0.2
        y_label_string = "Graph Density"
    elseif group_by == "cutoff"
        # bins = 1.0:-0.05:0.2
        bins = 0.5:-0.05:0.2
        # y_label_string = "Cutoff"
        y_label_string=" Normalized Cutoff Distance"
    end

    for (i,graph_size) in enumerate(graph_sizes)
        # df_sys_nc_pg = collect_results(datadir("exp_raw_e_new/$subresult_path/$graph_size"))
        df = collect_results("C:/age-inf-e/$(subresult_path)_07/$graph_size") 
        # df_pg = collect_results(datadir("exp_raw_e_new/$pg_path/$graph_size"))
        
        # df =  collect(df_pg, df_sys_nc) # 尝试合并 每个里面都有三个 

        gdfs = groupby(df, :cutoff, sort=true)  # 同一个cuttoff下

        g_x_values = []
        g_y_values = []
        g_z_values = []


        for gdf in gdfs # 同一个cuttoff下

            # sys_nc_pg_gdf = groupby(gdf, :max_inflight, sort=true)
            sys_inf_gdf = groupby(gdf, :max_inf, sort = true)

            # @info sys_inf_gdf
            sys_df = sys_inf_gdf[1]
            # @info sys_df
            nc_df = sys_inf_gdf[2]

            # @info graph_size
            # nc_pg_df  = groupby(nc_pg_df, :piggyback_enable, sort=true)
            # nc_df = nc_pg_df[2] # false true 因此第二个 是 pg


            ############## 下面就保持

            sort!(sys_df, [:graph_id])
            sort!(nc_df, [:graph_id])


            max_common_sample = min(nrow(sys_df),nrow(nc_df))
            # gains = [sys_df.aoi_components[i][1]/nc_df.aoi_components[i][1] for i = 1:max_common_sample]
            # gains = [sys_df.aoi_from_direct[i][1]/nc_df.aoi_from_dist[i] for i = 1:max_common_sample]
            gains = [sys_df.aoi_from_direct[i][1]/nc_df.aoi_from_direct[i][1] for i = 1:max_common_sample]
            # gains = [sys_df.aoi_from_dist[1]/nc_df.aoi_from_dist[i] for i = 1:max_common_sample]
            
            if group_by == "density"
                append!(g_y_values,sys_df.cutoff[1:max_common_sample])
            elseif group_by == "cutoff"
                append!(g_y_values, rounds.([sys_df[i,:].cutoff for i in 1:max_common_sample],0.05))
            end
            append!(g_z_values, (gains.-1).*100)
        end
        gp_z_values = []
        gp_sample_values = []
        for bin in bins
            samples_indices = findall(g_y_values .== bin)
            if !isempty(samples_indices)
                append!(gp_z_values, mean(g_z_values[samples_indices]))
                append!(gp_sample_values, length(g_z_values[samples_indices]))
            else
                append!(gp_z_values, 0)
                append!(gp_sample_values, 0)
            end

        end

        push!(x_values, graph_size)
        push!(y_values, bins)
        push!(z_values, gp_z_values)
        push!(sample_values, gp_sample_values)


    end
    
    x_dim = length(x_values)
    y_dim = maximum([length(z) for z in z_values])
    
    data = zeros(y_dim,x_dim)
    labels = zeros(y_dim,x_dim)
    for x in 1:length(z_values)
        for y in 1:length(z_values[x])
            data[y,x] = z_values[x][y]
            labels[y,x] = sample_values[x][y]
        end
    end
    fig, ax = PyPlot.subplots()
    xlabel("\$N\$", fontsize = 12)#,fontfamily  = "Times New Roman")
    ylabel(y_label_string, fontsize = 12)#,fontfamily  = "Times New Roman")
    im = ax.imshow(data,cmap="YlOrBr")
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.94, anchor=(0,0), pad=0.06)
    # cbar.ax.set_title("\$\\mathcal{G}_{inf}~[\\%]\$", fontsize = 10.5)#,fontfamily  = "Times New Roman")
    # cbar.ax.set_title("\$\\mathcal{G}~[\\%]\$", fontsize = 10.5)#,fontfamily  = "Times New Roman")
    cbar.ax.set_title("\$[\\%]\$", fontsize = 10.5)#,fontfamily  = "Times New Roman")

    ax.set_xticks((1:x_dim) .-1, labels=x_values)
    # label_set = string.(collect(1.0:-0.05:0.2))
    label_set = string.(collect(0.5:-0.05:0.2))
    label_set[2:2:end] .= ""
    ax.set_yticks((1:1:y_dim) .-1, labels=label_set, fontsize = 12)#,fontfamily  = "Times New Roman")

    #######################################################################################################
    for i in 1:x_dim
        for j in 1:y_dim
            # color = data[j,i]/maximum(data) > 0.75 ? [254,254,254]./255 : [0.0,0.0,0.0]
            color = data[j,i]/maximum(data) > 0.85 ? [254,254,254]./255 : [0.0,0.0,0.0]
            color = "#$(hex(RGB(color[1],color[2],color[3])))"
            ax.text(i-1, j-1, "$(trunc(Int8,data[j, i]))", ha="center", va="center", color=color, fontsize=8, linespacing=1.2)
        end
    end
    # im.set_clim(0,75)
    im.set_clim(0,85)
    # title("$title_string", fontsize = 12)#,fontfamily  = "Times New Roman")
    # title("\$\\mathcal{G}~[\\%]\$", fontsize = 12)#,fontfamily  = "Times New Roman")
    title("\$\\mathcal{G}\$", fontsize = 12)#,fontfamily  = "Times New Roman")
    savefig(plotsdir("inf-sys-sysinf-graphsize_degree_gain-$save_title_string-$group_by-big.png"),bbox_inches="tight")
    savefig(plotsdir("inf-sys-sysinf-graphsize_degree_gain-$save_title_string-$group_by-big.eps"),bbox_inches="tight")
end


pl_graphsize_degree_gain([12,14,16,18,20], "cutoff", "geometric_graphs")
# pl_graphsize_degree_gain([8,10,12], "cutoff", "geometric_graphs")
# pl_graphsize_degree_gain([8,10,12,14,16], "cutoff", "geometric_graphs")
# pl_graphsize_degree_gain([8,10], "density", "erdos_renyi")
# pl_graphsize_degree_gain([8,12,14,16,18,20], "density", "erdos_renyi")
# pl_graphsize_degree_gain([8,10,12,14,16,18,20], "density", "erdos_renyi")
# pl_graphsize_degree_gain([10], "density", "erdos_renyi")