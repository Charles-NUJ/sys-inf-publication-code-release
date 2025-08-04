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

using Graphs
using Combinatorics
# 检查一个候选的最小连通支配集（MCDS）是否覆盖了所有节点，并且该候选子图是连通的。
function is_connected_dominating_set(full_graph, mcds_nodes, mcds_candidate_graph)
    N = nv(full_graph) # 返回向量数
    #计算被覆盖的节点数，包括MCDS节点及其邻居
    nodes_covered = size(unique(vcat(mcds_nodes, reduce(vcat, [neighbors(full_graph, mcds_node) for mcds_node in mcds_nodes]))), 1)
    # 如果覆盖的节点数等于总节点数且候选图是连通的，则返回true
    if nodes_covered == N && is_connected(mcds_candidate_graph)
        return true
    else
        return false
    end
end
#将一个节点子集转换为子图，方法是从连接矩阵中提取这些节点的连接信息，生成新的子图。
function node_set_to_subgraph(connectivity, node_subset)
    node_subset = collect(node_subset) # 将节点子集转换为向量
    d = Array(zeros(size(node_subset, 1), size(node_subset, 1)))# 初始化一个零矩阵
    for (i_id, i) in enumerate(node_subset)
        for (j_id, j) in enumerate(node_subset)
            d[j_id, i_id] = connectivity[j, i]# 填充邻接矩阵
            d[i_id, j_id] = connectivity[j, i]
        end
    end
    return Graph(d)
end
#根据节点数量决定使用简化方法还是完整方法来获取节点类型。节点数少于等于16时使用完整方法，否则使用简化方法。
function get_node_types(connectivity)
    N = size(connectivity, 1)
    if N <= 16
        return get_node_types_mcds(connectivity) # determine leaves, pseudo-leaves (Farazi paper) and non-leaves, as well as MCDS size 确定 叶 伪叶 非叶 和 最小连通支配集的大小
    else
        return get_node_types_reduced(connectivity) # only finds the "true leaves" and consideres all other nodes non-leaves, does not calculate mcds
    end
end
# 简化版本，只找到真叶节点，并将其他节点视为非叶节点，不计算最小连通支配集（MCDS）。
function get_node_types_reduced(connectivity)
    non_leaves = []
    true_leaves = []
    pseudo_leaves = []
    mcds = []
    best_cds_size = 0
    graph = Graph(connectivity)
    N = nv(graph)
    nodes = 1:N
    for n in nodes
        if sum(connectivity[n, :]) == 1
            append!(true_leaves, n)# 如果节点度为1，则为真叶节点
        else
            append!(non_leaves, n)# 否则为非叶节点
        end
    end
    return best_cds_size, non_leaves, true_leaves, pseudo_leaves, mcds
end
#完整版本，找到叶、伪叶、非叶节点，并计算最小连通支配集（MCDS）的大小和组成。
function get_node_types_mcds(connectivity)
    non_leaves = [] # 非叶
    true_leaves = [] # 叶
    pseudo_leaves = [] # 伪叶

    graph = Graph(connectivity)
    # println(graph)
    N = nv(graph)
    nodes = 1:N

    mcds_candidates = combinations(nodes) # 所有元素的组合
    # println(mcds_candidates)
    best_cds_size = N
    mcds = []
    # Find a single minimum connected dominating set  # 找到一个最小连通支配集
    for mcds_candidate in mcds_candidates
        mcds_candidate = collect(mcds_candidate)
        mcds_candidate_graph = node_set_to_subgraph(connectivity, mcds_candidate)
        if is_connected_dominating_set(graph, mcds_candidate, mcds_candidate_graph)
            if length(mcds_candidate) < best_cds_size
                best_cds_size = length(mcds_candidate)
            end
        end
    end
    # Check all CDS, check if MCDS, union mcds nodes 联合mcds节点
    # # 检查所有CDS，找到MCDS，并将MCDS节点并入非叶节点集合
    for mcds_candidate in mcds_candidates
        mcds_candidate = collect(mcds_candidate)
        mcds_candidate_graph = node_set_to_subgraph(connectivity, mcds_candidate)
        if is_connected_dominating_set(graph, mcds_candidate, mcds_candidate_graph)
            if length(mcds_candidate) == best_cds_size
                union!(non_leaves, collect(mcds_candidate))
                mcds = deepcopy(mcds_candidate)
            end
        end
    end
    true_leaves = findall(vec(sum(connectivity, dims=1) .== 1)) #找到所有真叶节点的索引
    pseudo_leaves = setdiff(1:N, true_leaves, non_leaves) # # 找到伪叶节点（在N中但不在真叶和非叶集合中）
    return best_cds_size, non_leaves, true_leaves, pseudo_leaves, mcds
end
