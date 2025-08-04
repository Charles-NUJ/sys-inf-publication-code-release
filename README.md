
Any references to this code, protocol, or results shall be realized by citing the original paper from IEEE Xplore.

## Preface
Our work is based on the networking method of https://ieeexplore.ieee.org/document/10279069

The original code repository for the scientific publication is 'Optimizing Age of Information in Status Update Systems using Network Coding: A Graph Search Approach' by Fisser, Leonard and Timm-Giel, Andreas published at the IEEE International Conference on Communications 2023 held in Rome, Italy.

The modified file is mainly in state_utilities_src_inf_v3.jl and incremental_a_star_inf_v3.jl

After the paper is accepted, we will disclose all the project codes.

## Getting Started
1. Download and install [JuliaLang](https://julialang.org/downloads/oldreleases/) Version 1.7.3.
2. Clone the repository to your local computer.
3. Start JuliaLang and navigate to the project's root directory.
4. Install project dependencies using:
```julia
pkg> activate .
pkg> instantiate
```

```julia
# 动作的数据结构
single_type_and_action = Dict(
    "action_type" => "single",
    "action" => single_action,
)
sim_type_and_action = Dict(
    "action_type" => "sim",
    "action" => sim_action,
    # "sim_tx_number" => length(sim_action),
)
# 负载的数据结构
payload_dict = Dict(
    "payload_type" => "single",
    "tx" => action[1],
    "payload" => payload,
)
payload_dict = Dict(
    "payload_type" => "sim",
    "tx_list" => tx_list,
    "payload_sim" => payload_s,
)

# 进一步将Key简化为
# "action_type" "at"
# "action" "a"
# "payload_type" "pt"
# "payload" ""payload_sim" "p"
# "tx” "tx_list"  "t"

# 邻居的数据结构 在利用邻居信息的时候可能用到
UpdateNeibors_all_simu_dict_temp = Dict(
    "neighbor_type" => "single",
    "UpdateNeibors" => UpdateNeibors,
)
UpdateNeibors_all_simu_dict_temp = Dict(
    "neighbor_type" => "sim",
    "UpdateNeibors" => UpdateNeibors_all_simu,
)
```



 [fig.four_line.pdf](git_figure\fig.four_line.pdf) 

![inf_aoi_from_direct_C_e](git_figure\inf_aoi_from_direct_C_e.png)



![aoi-CMP_sys_inf_nc_inf_pool2A](git_figure\aoi-CMP_sys_inf_nc_inf_pool2A.png)

![inf-sys-sysinf-graphsize_degree_gain-erod_renyi_graphs-density-big](git_figure\inf-sys-sysinf-graphsize_degree_gain-erod_renyi_graphs-density-big.png)

![inf-sys-sysinf-graphsize_degree_gain-erod_renyi_graphs-density-big-T](git_figure\inf-sys-sysinf-graphsize_degree_gain-erod_renyi_graphs-density-big-T.png)

![inf-sys-sysinf-graphsize_degree_gain-geometric_graphs-cutoff-big](git_figure\inf-sys-sysinf-graphsize_degree_gain-geometric_graphs-cutoff-big.png)

![inf-sys-sysinf-graphsize_degree_gain-geometric_graphs-cutoff-big-T](git_figure/inf-sys-sysinf-graphsize_degree_gain-geometric_graphs-cutoff-big-T.png)
