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

using DataFrames, Statistics, PyPlot, Graphs, GraphIO

include("$(srcdir())/utilities/topology_utilities.jl")
include("$(srcdir())/utilities/aoi_utilities.jl")




function min_avg_min_id(pool)
    # d4 = pool[1]
    # d5 = pool[2]
    # d6 = pool[3]
    # d7 = pool[3]
    
    # return [minimum(a),mean(t),maximum(d),max_index[2]]

    for d in pool
        max_index = findmax(d)
        min_index = findmin(d)
        @info "$(minimum(d)*100), $min_index, $(mean(d)*100), $(maximum(d)*100), $max_index"
    end
end


function pl_aoi_components(graph_sizes)
    
    gain_update = "#CCBB44"
    gain_encoding = "#EE6677"
    gain_full = "#4477AA"

    full_set_gains = []
    full_set_t_diff = []

    full_set_gains_bigger = []
    full_set_gains_all = []
    full_set_gains_bigger_all = []


    T_set_gains_all = []
    T_set_gains_all_bigger_all = []

    D_set_gains_all = []
    D_set_gains_all_bigger_all = []

    a_pool=[]
    t_pool=[]
    d_pool=[]

    # fig, ax = PyPlot.subplots(1,4) # 四幅图
    fig, ax = PyPlot.subplots(2,2)

    fig.tight_layout(rect=(0.045,0,1,0.95))
    fig.subplots_adjust(wspace=0.05, hspace=0.2)
    
    # ax = [ax[1],ax[3],ax[2],ax[4]]
    # ax = [ax[1],ax[2],ax[3],ax[4]]
    ax = [ax[1],ax[3],ax[2],ax[4]]

    plot_labeling = ["a)","b)","c)","d)"]

    ## Set up LaTeX fonts
    plt.rcParams["text.usetex"] = true
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serifx"] =  ["Times New Roman"]
    plt.rcParams["font.size"] = 30
    # plt.rcParams['pdf.fonttype'] = 42
    plt.rc("hatch", linewidth=3.25)
    Fontsize=14
    for (i,graph_size) in enumerate(graph_sizes)

        # @info i graph_size
        # df = collect_results(datadir("exp_raw/$graph_size"))

        # df = collect_results("D:/AoICode/ICC2023/age-optimal-multisource-flooding/data/exp_raw_sys_inf/pool_sum+dist_both/$graph_size")
        # df = collect_results("C:/age-optimal-multisource-flooding/data/exp_raw_sys_inf/pool/$graph_size")
        # df = collect_results("C:/age-optimal-multisource-flooding/data/exp_raw_sys_inf_e/pool_3/$graph_size")
        df = collect_results("C:/age-optimal-multisource-flooding/data/exp_raw_sys_inf_e/pool_good/$graph_size")
        
        # df = collect_results("D:/AoICode/ICC2023/age-optimal-multisource-flooding/data/exp_raw_sys_inf/pool_sum+dist/$graph_size")
        # df = collect_results("D:/AoICode/ICC2023/age-optimal-multisource-flooding/data/exp_raw_sys_inf/pool_sum+dist_sys_sys_inf_no/$graph_size")



        
        sys_inf_gdf = groupby(df, :max_inf, sort=true)  # 按照这个一排， 分为两组 max_inflight 为1和graph_sizes(节点数n) 4 5 6 7
 
        sys_df = sys_inf_gdf[1] # 编码度为1
        inf_df = sys_inf_gdf[2] # 编码最大度为n
 



        sort!(sys_df, [:graph_id]) # 按照 n中有多少个子图排序 # 4有     共6拓扑    2是线性拓扑 # 5节点   共21拓扑   9是线性拓扑 # 6节点   共112拓扑  # 7节点   共853拓扑
        sort!(inf_df, [:graph_id]) 
        
        if graph_size==4
            @info sys_df
            @info inf_df
         end


        g_aoi_full = [sys[1] for sys in sys_df.aoi_from_direct]./[inf[1] for inf in inf_df.aoi_from_direct] .- 1
        g_aoi_update = ([sys[2] for sys in sys_df.aoi_from_direct] - [inf[2] for inf in inf_df.aoi_from_direct]) ./ [inf[1] for inf in inf_df.aoi_from_direct]
        g_aoi_encoding = ([sys[3] for sys in sys_df.aoi_from_direct] - [inf[3] for inf in inf_df.aoi_from_direct]) ./ [inf[1] for inf in inf_df.aoi_from_direct]
    
        # 对sys 采用 aoi_components 或者 aoi_from_direct
        # 对sys-inf 采用 aoi_from_dist

        g_aoi_full_bigger = []

        a_tmp=[]
        t_tmp=[]
        d_tmp=[]


        for i in 1:length(g_aoi_full)
            if g_aoi_full[i] > 0
                 push!(g_aoi_full_bigger, g_aoi_full[i] )
                 push!(full_set_gains_bigger_all,g_aoi_full[i]) # 保存所有的比较大的gains_bigger
                 push!(T_set_gains_all_bigger_all,g_aoi_update[i] ) # 保存所有 增益
                 push!(D_set_gains_all_bigger_all,g_aoi_encoding[i] ) # 保存所有 增益

                 if graph_size==4
                    @info "index = $i"
                    @info  sys_df.schedule[i]
                    @info  inf_df.schedule[i]
                 end
            end
            push!(full_set_gains_all,g_aoi_full[i] ) # 保存所有 增益
            push!(T_set_gains_all,g_aoi_update[i] ) # 保存所有 增益
            push!(D_set_gains_all,g_aoi_encoding[i] ) # 保存所有 增益

            push!(a_tmp,g_aoi_full[i] ) # 保存所有 增益
            push!(t_tmp,g_aoi_update[i] ) # 保存所有 增益
            push!(d_tmp,g_aoi_encoding[i] ) # 保存所有 增益
        end
  

        push!(a_pool,a_tmp) # 保存所有 增益
        push!(t_pool,t_tmp) # 保存所有 增益
        push!(d_pool,d_tmp) # 保存所有 增益


        
        sorter = sortperm(g_aoi_full) # 排序，数据中要是有两个相同的，按照其出现的先后顺序放置
        g_aoi_full = g_aoi_full[sorter]
        # g_aoi_full_capped = g_aoi_full_capped[sorter]
        g_aoi_update = g_aoi_update[sorter]
        g_aoi_encoding = g_aoi_encoding[sorter]

        # g_aoi_full = g_aoi_full_bigger
        # sorter = sortperm(g_aoi_full) # 排序，数据中要是有两个相同的，按照其出现的先后顺序放置
        # g_aoi_full = g_aoi_full[sorter]

        g_bar_positions = collect(1:length(g_aoi_full))
        if i > 2
            bar_width = 1.0
        else
            bar_width = 0.8
        end
        ratio = 0.66
        ax[i].bar(g_bar_positions, g_aoi_update.*100, color=gain_update, width=bar_width, edgecolor = gain_update)        
        
        if i > 2
            ratio = 1.0
            ax[i].bar(g_bar_positions.+(bar_width*(1-ratio))/2, g_aoi_full.*100 , color=gain_full, width=bar_width*ratio, edgecolor = gain_full)
            ax[i].bar(g_bar_positions, g_aoi_encoding.*100, color=gain_encoding, width=bar_width, edgecolor = gain_encoding)
        else
            ax[i].bar(g_bar_positions.+(bar_width*(1-ratio))/2, g_aoi_full.*100 , color=gain_full, width=bar_width*ratio, edgecolor = gain_full)
            ax[i].bar(g_bar_positions, g_aoi_encoding.*100, color=gain_encoding, width=bar_width, edgecolor = gain_encoding)
        end

        plt.sca(ax[i])
        xticks([])

        if i==1 || i==2
        ylim([-10,20])
        end

        if i==3 || i==4
            ylim([-10,75])
        end

        xlabel("$(plot_labeling[i]) \$N=$graph_size\$",fontsize=12)
        if mod(i,2) == 1
            ylabel("Gain [%]")
        else
            # plt.setp(ax[i], yticklabels=[])
            plt.setp(ax[i], yticklabels=[])
        end
        ax[i].set_axisbelow(true)
        ax[i].grid(which="major",  color="#e7e7e7ff", linestyle="-", linewidth=0.5)
        
        if i==1 || i==2
            yticks(collect(-10:5:20), fontsize=12)
          end
        if i==3 || i==4
          yticks(collect(-10:10:75), fontsize=12)
        end

        append!(full_set_gains,g_aoi_full)
        # append!(full_set_t_diff,inf_df.t-inf_df.t_star)
        append!(full_set_t_diff,inf_df.t-inf_df.t_star)
        append!(full_set_gains_bigger,g_aoi_full_bigger)
    end
    
    lines = []
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_full))
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_update))
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_encoding))
    
    # labels = ["\$\\mathcal{G}\$"]

    labels = [
        "\$\\mathcal{G}\$",
            "\$\\mathcal{T}\$",
              "\$\\mathcal{W}+\\mathcal{L}\$",
            #   "\$\\mathcal{G}\$"
            #   "\$U^{2}_{\\Delta}+C^{2}_{\\Delta}\$"
              ]
    
    fig.legend(lines, labels, loc = (0.28, 0.925), ncol=4, edgecolor="black", fontsize=12)

    savefig(plotsdir("inf_aoi_from_direct_C_e.png"),bbox_inches="tight")
    savefig(plotsdir("inf_aoi_from_direct_C_e.eps"),bbox_inches="tight")
    # savefig(plotsdir("inf_aoi_from_direct.pdf"),bbox_inches="tight")
    # savefig(plotsdir("inf_aoi_from_direct.svg"),bbox_inches="tight")
    println("Average gain: $(mean(full_set_gains)*100) %")
    println("Maximum gain: $(maximum(full_set_gains)*100) %")
    # println("In $(count(full_set_t_diff.>0)) out of $(size(full_set_t_diff)[1]) schedules, non t_star was optimal.")
    println("Average gain in bigger : $(mean(full_set_gains_bigger)*100) %")


    println("全部有个数据: $(length(full_set_gains_all))")
    # println("G全部的增益平均为: $(mean(full_set_gains_all)*100) %")
    # println("G为正的增益的平均为: $(mean(full_set_gains_bigger_all)*100) %")
    println("G增益为正的有:  $(length(full_set_gains_bigger_all))")
    println("G全部的最小为$(minimum(full_set_gains_all)*100) %，最大为$(maximum(full_set_gains_all)*100) %，平均为: $(mean(full_set_gains_all)*100) %")
    println("T全部的最小为$(minimum(T_set_gains_all)*100) %，最大为$(maximum(T_set_gains_all)*100) %，平均为: $(mean(T_set_gains_all)*100) %")
    println("D全部的最小为$(minimum(D_set_gains_all)*100) %，最大为$(maximum(D_set_gains_all)*100) %，平均为: $(mean(D_set_gains_all)*100) %")



    len = length(full_set_gains_all)
    println("T全部变小的总数比例为$((length(findall(T_set_gains_all.>0))/len)*100) %，不变为$((length(findall(T_set_gains_all.==0))/len)*100) %，变大为: $((length(findall(T_set_gains_all.<0))/len)*100) %")
    println("D全部变小的总数比例为$((length(findall(D_set_gains_all.>0))/len )*100) %，不变为$((length(findall(D_set_gains_all.==0))/len)*100) %，变大为: $((length(findall(D_set_gains_all.<0))/len)*100) %")
    println("A全部变小的总数比例为$((length(findall(full_set_gains_all.>0))/len)*100) %，不变为$((length(findall(full_set_gains_all.==0))/len)*100) %，变大为: $((length(findall(full_set_gains_all.<0))/len)*100) %")

    println("T全部变小的总数为$((length(findall(T_set_gains_all.>0))/1)*1) ，不变为$((length(findall(T_set_gains_all.==0))/1)*1) ，变大为: $((length(findall(T_set_gains_all.<0))/1)) ")
    println("D全部变小的总数为$((length(findall(D_set_gains_all.>0))/1 )*1) ，不变为$((length(findall(D_set_gains_all.==0))/1)*1) ，变大为: $((length(findall(D_set_gains_all.<0))/1)) ")
    println("A全部变小的总数为$((length(findall(full_set_gains_all.>0))/1)*1) ，不变为$((length(findall(full_set_gains_all.==0))/1)*1) ，变大为: $((length(findall(full_set_gains_all.<0))/1)) ")
    # println("D增益为正的有:  $(length(full_set_gains_bigger_all))")

    @info "_____________________________________________________________________________________________________"

    println("G为正的比例为 $(( length(full_set_gains_bigger_all) /len )*100) % ")

    println("G为正中G最小为$(minimum(full_set_gains_bigger_all)*100) %，最大为$(maximum(full_set_gains_bigger_all)*100) %，平均为: $(mean(full_set_gains_bigger_all)*100) %")
    println("G为正中T最小为$(minimum(T_set_gains_all_bigger_all)*100) %，最大为$(maximum(T_set_gains_all_bigger_all)*100) %，平均为: $(mean(T_set_gains_all_bigger_all)*100) %")
    println("G为正中D最小为$(minimum(D_set_gains_all_bigger_all)*100) %，最大为$(maximum(D_set_gains_all_bigger_all)*100) %，平均为: $(mean(D_set_gains_all_bigger_all)*100) %")

    len = length(full_set_gains_bigger_all)
    println("G为正中T变小为$((length(findall(T_set_gains_all_bigger_all.>0))/len )*100) %，不变为$((length(findall(T_set_gains_all_bigger_all.==0))/len)*100) %，变大为: $((length(findall(T_set_gains_all_bigger_all.<0))/len)*100) %")
    println("G为正中D变小为$((length(findall(D_set_gains_all_bigger_all.>0))/len )*100) %，不变为$((length(findall(D_set_gains_all_bigger_all.==0))/len)*100) %，变大为: $((length(findall(D_set_gains_all_bigger_all.<0))/len)*100) %")

    println("G为正中T 变小的总数为$((length(findall(T_set_gains_all_bigger_all.>0))/1)*1) ，不变为$((length(findall(T_set_gains_all_bigger_all.==0))/1)*1) ，变大为: $((length(findall(T_set_gains_all_bigger_all.<0))/1)) ")
    println("G为正中D 变小的总数为$((length(findall(D_set_gains_all_bigger_all.>0))/1 )*1) ，不变为$((length(findall(D_set_gains_all_bigger_all.==0))/1)*1) ，变大为: $((length(findall(D_set_gains_all_bigger_all.<0))/1)) ")
    println("G为正中A 变小的总数为$((length(findall(full_set_gains_bigger_all.>0))/1)*1) ，不变为$((length(findall(full_set_gains_bigger_all.==0))/1)*1) ，变大为: $((length(findall(full_set_gains_bigger_all.<0))/1)) ")




    @info "对于每组数据"
    @info "age" 
    min_avg_min_id(a_pool)

    @info "t" 
    min_avg_min_id(t_pool)

    @info "d"  
    min_avg_min_id(d_pool)
end

graph_sizes = [4,5,6,7] # 
pl_aoi_components(graph_sizes)

# df = collect_results("D:/AoICode/ICC2023/age-optimal-multisource-flooding/data/exp_raw_sys_inf/pool/3")

# @info df

# The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.
# Average gain: 3.7053081542918145 %
# Maximum gain: 72.15447154471546 %
# Average gain in bigger : 6.832092358842903 %
# 全部有个数据: 992
# G增益为正的有:  538
# G全部的最小为0.0 %，最大为72.15447154471546 %，平均为: 3.7053081542918163 %
# T全部的最小为0.0 %，最大为72.56097560975611 %，平均为: 4.341385487016322 %
# D全部的最小为-5.633802816901412 %，最大为0.4040404040404045 %，平均为: -0.6360773327245082 %
# T全部变小的总数比例为54.233870967741936 %，不变为45.766129032258064 %，变大为: 0.0 %
# D全部变小的总数比例为0.3024193548387097 %，不变为59.2741935483871 %，变大为: 40.42338709677419 %
# A全部变小的总数比例为54.233870967741936 %，不变为45.766129032258064 %，变大为: 0.0 %
# [ Info: _____________________________________________________________________________________________________
# G为正的比例为 54.233870967741936 %
# G为正中G最小为1.3513513513513375 %，最大为72.15447154471546 %，平均为: 6.832092358842903 %
# G为正中T最小为3.5234899328859064 %，最大为72.56097560975611 %，平均为: 8.004933834795896 %
# G为正中D最小为-5.633802816901412 %，最大为0.4040404040404045 %，平均为: -1.1728414759529964 %
# G为正中T变小为100.0 %，不变为0.0 %，变大为: 0.0 %
# G为正中D变小为0.5576208178438662 %，不变为24.907063197026023 %，变大为: 74.53531598513011 %
# [ Info: 对于每组数据
# [ Info: age
# [ Info: 0.0, (0.0, 1), 0.42735042735042955, 2.564102564102577, (0.025641025641025772, 2)
# [ Info: 0.0, (0.0, 1), 1.211927146242934, 7.106598984771595, (0.07106598984771595, 9)
# [ Info: 0.0, (0.0, 1), 3.1024455225581526, 63.84615384615384, (0.6384615384615384, 47)
# [ Info: 0.0, (0.0, 1), 3.8689064688109793, 72.15447154471546, (0.7215447154471546, 422)
# [ Info: t
# [ Info: 0.0, (0.0, 1), 1.2820512820512822, 7.6923076923076925, (0.07692307692307693, 2)
# [ Info: 0.0, (0.0, 1), 2.0388310356772315, 10.152284263959391, (0.10152284263959391, 9)
# [ Info: 0.0, (0.0, 1), 3.9754364118596763, 63.46153846153847, (0.6346153846153847, 47)
# [ Info: 0.0, (0.0, 1), 4.467640991266566, 72.56097560975611, (0.7256097560975611, 422)
# [ Info: d
# [ Info: -5.128205128205128, (-0.05128205128205127, 2), -0.8547008547008547, 0.0, (0.0, 1)
# [ Info: -3.164556962025316, (-0.03164556962025316, 4), -0.8269038894342964, 0.0, (0.0, 1)
# [ Info: -5.633802816901412, (-0.05633802816901412, 4), -0.8729908893015244, 0.3846153846153833, (0.003846153846153833, 47)
# [ Info: -3.9215686274509785, (-0.03921568627450978, 4), -0.5987345224555868, 0.4040404040404045, (0.004040404040404045, 592)