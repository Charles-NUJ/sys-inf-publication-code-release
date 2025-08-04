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

using Graphs, Combinatorics
using JuMP, GLPK, LinearAlgebra#, Gurobi
#using Gurobi # Consider using Gurobi for large problems
include("$(srcdir())/utilities/topology_utilities.jl")


# 计算网络编码洪范所需要的最小节点容量
function min_schedule_length_nc(connectivity, inactive_nodes = [])
    """Returns the required node capacities for full network coded flooding. Taken from Jörg Widmer paper."""
    model = Model(GLPK.Optimizer) # 创建优化模型
    g = Graph(connectivity) # 创建图对象
    N = nv(g)
    @variable(model, C[1:N], Int) # 定义节点容量变量 C
    for c in C
        @constraint(model, c >= 0) #设置约束 C>=0
    end
    for cut in combinations(collect(1:N)) # 遍历所有可能的节点切割
        S = cut
        S_prime = setdiff(collect(1:N),S)# 计算切割后的节点集合S 和 S_prime
        if size(S,1) < N 
            cut_mask = zeros(Int,N)
            cut_mask[S] .= 1
            cut_mask[S_prime] .= 2
            Ns = unique(cat([[n.src n.dst] for n in karger_cut_edges(g, cut_mask)]...,dims=2))# 计算切割边的节点
            intersect!(Ns,S)
            setdiff!(S,inactive_nodes)
            @constraint(model,sum(C[Ns]) >= size(S,1)) # 添加约束，确保总容量 》=切割集合的大小
        end
    end
    @objective(model, Min, sum(C)) # 设置目标函数：最小化节点容量总和
    optimize!(model)
    return JuMP.value.(C) # 返回最优节点容量
end


# 建立一个动态优化问题，可以调整活跃节点的集合
function get_dynamic_opt_problem(connectivity)
    """Returns a dynamic optimization problem, for which the set of active nodes can be adjusted."""
    model = Model(GLPK.Optimizer)
    #model = Model(Gurobi.Optimizer) #For large problems, Gurobi's performance is required
    set_silent(model)
    g = Graph(connectivity)
    N = nv(g)
    model[:C] = @variable(model, C[1:N], Int) # 容量变量 C
    model[:node_active] = @variable(model, node_active[1:N], Int)# 活跃状态变量 node_active
    for c in C
        @constraint(model, c >= 0) # 节点容量的非负约束
    end
    for cut in combinations(collect(1:N)) # 便利所以可能的节点切割
        S = cut
        S_prime = setdiff(collect(1:N),S) # 计算切割后的节点集合S 和 S_prime
        if size(S,1) < N 
            cut_mask = zeros(Int,N)
            cut_mask[S] .= 1
            cut_mask[S_prime] .= 2
            Ns = unique(cat([[n.src n.dst] for n in karger_cut_edges(g, cut_mask)]...,dims=2))  # 计算切割边的节点
            intersect!(Ns,S)
            @constraint(model,sum(C[Ns]) >= sum(node_active[S])) # 添加约束，确保切割集合中的活跃节点总数 >= 该切割上的总容量
        end
    end
    @objective(model, Min, sum(C))  # 设置目标函数：最小化节点容量总和
    return model  # 返回构建好的优化模型
end