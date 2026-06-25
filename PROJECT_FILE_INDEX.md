# 项目代码与数据文件说明

本文档按研究流程整理项目中的主要代码、数据和绘图文件。它用于快速定位文件职责，不替代源码注释和论文中的数学定义。除特别说明外，Julia 脚本均建议在项目根目录运行：

```powershell
julia --project=. path\to\script.jl
```

## 1. 运行环境与根目录文件

- `Project.toml`：Julia 项目依赖声明。
- `Manifest.toml`：Julia 依赖锁定版本。
- `README.md`：原始项目说明。
- `LICENSE`：开源许可。
- `PROJECT_FILE_INDEX.md`：当前文件，项目代码与数据索引。
- `log1 无编码.txt`：全局干扰、无编码调度日志，用于理解动作和状态输出格式。
- `checkresult.jl`、`checkresult_leavf.jl`：历史结果检查脚本。
- `tmp_pan5_compact.jl`：五节点 pan 小例子的临时检查脚本。

## 2. 核心源码目录

### 2.1 基础图搜索与状态工具

- `src/solvers/a_star_src.jl`：全局干扰模型下的精确 A* 搜索；每个时隙只允许一个节点广播。
- `src/solvers/incremental_a_star_src.jl`：全局干扰模型下的增量 A* / prefix-and-rollback 搜索。
- `src/solvers/nc_capacities.jl`：网络编码相关容量和辅助计算，用于和文献结果对比，不是本文大规模主算法的核心。
- `src/utilities/state_utilities_src.jl`：全局干扰模型的状态、动作、状态转移工具。
- `src/utilities/aoi_utilities_src.jl`：全局干扰模型的 AoI、首次发送、首次接收、调度评价工具。
- `src/utilities/bin_mat_rref.jl`：二进制矩阵行最简形工具，主要服务网络编码状态判断。
- `src/utilities/topology_utilities.jl`：拓扑生成、读取和基础图工具。
- `src/utilities/quick_plot_graph.py`：快速绘制拓扑图的 Python 工具。

### 2.2 协议干扰模型源码

- `src/solvers/a_star_src_inf_v3.jl`：早期协议干扰 A* 实现，保留原始状态语义和输出格式。
- `src/solvers/a_star_src_inf_v8.jl`：协议干扰精确 A* 遍历版本。默认 `max_inflight=size(connectivity,1)`，可选保存小图遍历 trace。
- `src/solvers/a_star_src_inf_v9.jl`：从 v8 分离出的 heuristic/repair 版本，用于快速修复或获得可行解；函数名中仍保留部分历史 `v8` 命名。
- `src/solvers/incremental_a_star_inf_v3.jl`：协议干扰增量 A* 早期版本。
- `src/solvers/incremental_a_star_inf_v8.jl`：协议干扰增量 A* 的修正版本，用于排查 v3 的 pan-5 不一致问题。
- `src/solvers/incremental_a_star_inf_v9.jl`：协议干扰增量 A* 的可变参数版本，便于设置并发上限和 in-flight 上限。
- `src/utilities/state_utilities_src_inf_v3.jl`：协议干扰模型的状态、动作、并发可行性和状态转移工具。
- `src/utilities/aoi_utilities_src_inf_v3.jl`：协议干扰模型的 AoI 计算、首次接收矩阵和调度评价工具。

## 3. 小规模精确搜索与验证脚本

### 3.1 小图 A* 与 pan-5 验证

- `scripts/evaluations/five_node_pan_exact.jl`：五节点 pan 拓扑精确调度和 AoI 结果生成。
- `scripts/evaluations/validate_pan5_inf_astar.jl`：验证协议干扰 A* 在 pan-5 上的基本结果。
- `scripts/evaluations/validate_pan5_inf_astar_v8.jl`：验证 `a_star_src_inf_v8.jl` 的 pan-5 精确结果。
- `scripts/evaluations/validate_pan5_inf_astar_epsilon.jl`：检查 pan-5 A* 中 epsilon/容差设置对结果的影响。
- `scripts/evaluations/validate_pan5_astar_models.jl`：对比全局干扰 A* 与协议干扰 A* 的状态语义和结果。
- `scripts/evaluations/validate_pan5_age_first.jl`：在 pan-5 上运行 Age-first 多启动验证。
- `scripts/evaluations/validate_pan5_age_schedule_in_astar_state.jl`：把 Age-first 调度回放到 A* 状态机中，检查状态语义一致性。
- `scripts/evaluations/validate_pan5_incremental_inf.jl`：验证 `incremental_a_star_inf_v3.jl` 的 pan-5 结果。
- `scripts/evaluations/validate_pan5_incremental_v8_v9.jl`：比较增量 A* v8/v9 的 pan-5 结果。
- `scripts/evaluations/validate_pan5_top_e_incremental.jl`：通过 `top_e` 入口检查增量 A* 的 pan-5 运行流程。
- `scripts/evaluations/check_pan_node_types.jl`：检查 pan 拓扑节点类型和编号约定。

### 3.2 精确 MILP / BPC 验证

- `scripts/evaluations/small_exact_validation.jl`：小规模紧凑 MILP 与有限列 Dantzig-Wolfe 主问题等价性验证。
- `scripts/evaluations/small_branch_price_cut_validation.jl`：小规模 branch-price-and-cut 验证，包含动态定价、割、分支和外层周期搜索。
- `scripts/evaluations/pan5_fixed_bpc_validation.jl`：pan-5 固定周期 BPC / 有限列验证。
- `scripts/evaluations/uncoded_column_generation.jl`：无编码列生成建模与求解原型。
- `scripts/evaluations/uncoded_exact_100_star.jl`：无编码精确模型的历史检查脚本。

## 4. 中等规模增量搜索脚本

- `scripts/evaluations/top_e_src_inf_analysis_runner_inf.jl`：协议干扰增量 A* 的实验入口。
- `scripts/evaluations/top_e_src_inf_graph_analysis_inf.jl`：增量 A* 的拓扑读取、参数设置、结果保存逻辑。
- `scripts/evaluations/incremental_interference_astar_refactored.jl`：重构版增量协议干扰 A* 实验脚本。
- `scripts/evaluations/legacy_scalability_validation.jl`：将新的大规模启发式结果与历史中小规模增量 A* 结果进行对比。
- `data/exp_res/incremental_astar_20_results.txt`：20 节点增量 A* 的历史运行摘要。
- `data/exp_res/seven_node_scalability_validation.csv`：7 节点拓扑上的交叉验证结果。
- `data/exp_res/sixteen_node_euc_scalability_validation.csv`：16 节点欧氏图交叉验证结果。
- `data/exp_res/sixteen_node_renyi_scalability_validation.csv`：16 节点 Erdős-Rényi 图交叉验证结果。

## 5. 大规模集中式启发式调度

### 5.1 主算法入口

- `scripts/evaluations/interference_all_to_all_demo.jl`：协议干扰全对全调度核心构造器，包含图生成、冲突图、可行并发、AoI 评价等主流程。
- `scripts/evaluations/interference_algorithm_comparison.jl`：Coverage-first、Age-first、Dual-guided 三种集中式策略对比。
- `scripts/evaluations/age_first_large_scale.jl`：100/400/700/1000 等大规模 Age-first 运行与结果保存脚本。
- `scripts/evaluations/large_policy_stress_experiment.jl`：多策略大规模压力实验。
- `scripts/evaluations/scalability_statistics_ablation.jl`：固定半径、多随机种子、参数消融与可扩展性统计。
- `scripts/evaluations/validate_large_age_first.jl`：大规模 Age-first 调度合法性和结果检查。
- `scripts/evaluations/lower_bound_comparison_100.jl`：100 节点拓扑下多层下界、计算时间和松弛程度比较。
- `scripts/evaluations/systematic_dfs_baseline.jl`：复现文献中的 global-interference systematic DFS 基线。

### 5.2 大规模输出数据

- `data/exp_res/age_first_large_iterations.csv`：大规模 Age-first 每轮迭代结果。
- `data/exp_res/age_first_large_summary.csv`：大规模 Age-first 汇总结果。
- `data/exp_res/age_first_100_1000_summary.csv`：100/400/700/1000 节点摘要表。
- `data/exp_res/scalability_statistics_raw.csv`：固定半径多种子可扩展性原始记录。
- `data/exp_res/scalability_statistics_summary.csv`：可扩展性统计均值、方差、成功率、时间和内存。
- `data/exp_res/heuristic_ablation_raw.csv`、`data/exp_res/heuristic_ablation_raw_ablation.csv`：启发式参数消融原始记录。
- `data/exp_res/heuristic_ablation_raw_summary.csv`、`data/exp_res/heuristic_ablation_raw_ablation_summary.csv`：消融实验汇总。
- `data/exp_res/lower_bound_comparison_100.csv`：100 节点下界数值和运行时间。
- `data/exp_res/systematic_dfs_results.csv`：systematic DFS 基线结果。
- `data/exp_res/systematic_dfs_100_schedule.csv`：100 节点 systematic 调度序列。
- `data/exp_res/large_age_first_validation.csv`：大规模 Age-first 约束验证结果。
- `data/exp_res/large_age_first_validation_schedule.csv`：对应调度序列。

## 6. 大规模分布式部署仿真

- `scripts/evaluations/distributed_deployment_experiment.jl`：100 节点 Distributed-local 部署仿真，使用两跳 request/winner 竞争和本地估计。
- `scripts/evaluations/distributed_scalability_experiment.jl`：100/400/700/1000 节点 Distributed-local 可扩展性压力实验。
- `data/exp_res/distributed_deployment_100.csv`：100 节点分布式仿真逐轮结果。
- `data/exp_res/distributed_deployment_100_summary.csv`：100 节点分布式仿真最优摘要。
- `data/exp_res/distributed_deployment_100_schedule.csv`：100 节点分布式调度序列。
- `data/exp_res/distributed_scalability_raw.csv`：大规模分布式多尺寸原始结果。
- `data/exp_res/distributed_scalability_summary.csv`：大规模分布式汇总结果。
- `data/exp_res/distributed_scalability_best_N100_schedule.csv`：100 节点分布式最优调度序列。
- `data/exp_res/distributed_scalability_best_N400_schedule.csv`：400 节点分布式最优调度序列。
- `data/exp_res/distributed_scalability_best_N700_schedule.csv`：700 节点分布式最优调度序列。
- `data/exp_res/distributed_scalability_best_N1000_schedule.csv`：1000 节点分布式最优调度序列。

## 7. 服务器批处理脚本

这些脚本用于在服务器上批量读取 `C:\AoI_Inf\exp_res` 或本地 `data/exp_res` 中的 `.g6` 拓扑，并运行协议干扰 A* / 增量 A*。

- `scripts/evaluations/server_inf_batch_common.jl`：服务器批处理公共函数，包含环境变量读取、外部拓扑同步、参数打印等。
- `scripts/evaluations/server_top_src_inf_batch.jl`：批量运行协议干扰精确 A*。
- `scripts/evaluations/server_top_e_src_inf_batch.jl`：批量运行协议干扰增量 A*。
- `scripts/evaluations/top_src_inf_analysis_runner_inf_v8.jl`：协议干扰精确 A* v8 单实验入口。
- `scripts/evaluations/top_src_inf_graph_analysis_inf_v8.jl`：协议干扰精确 A* v8 的图读取和结果保存逻辑。
- `data/exp_res/graph3c.g6` 至 `data/exp_res/graph7c.g6`：小图 graph6 拓扑集合。
- `data/exp_res/graph_source.txt`：拓扑来源说明。

常用环境变量由 `server_inf_batch_common.jl` 统一解析，包括外部拓扑目录、图规模、图类型、图编号范围、是否覆盖已有结果等。

## 8. 约束验证与调度回放

- `scripts/evaluations/validate_saved_schedule_constraints.jl`：检查已保存调度是否满足协议干扰、因果传播、首次接收和 AoI 计算约束。
- `data/exp_res/saved_schedule_constraint_validation_100.csv`：100 节点调度约束验证摘要。
- `data/exp_res/saved_schedule_constraint_validation_best_100node_schedule.csv`：论文最优 100 节点调度回放验证。
- `data/exp_res/saved_schedule_constraint_validation_age_first_100_schedule.csv`：Age-first 100 节点调度回放验证。
- `data/exp_res/saved_schedule_constraint_validation_coverage_first_100_schedule.csv`：Coverage-first 100 节点调度回放验证。
- `data/exp_res/saved_schedule_constraint_validation_dual_guided_100_schedule.csv`：Dual-guided 100 节点调度回放验证。
- `data/exp_res/saved_schedule_constraint_validation_distributed_deployment_100_schedule.csv`：Distributed-local 100 节点调度回放验证。

## 9. 绘图与论文图表生成

### 9.1 论文 Fig. 1、Fig. 2、Fig. 8 和大规模图

- `scripts/figures/generate_figure1_pan.jl`：生成五节点 pan 拓扑对比图。
- `scripts/figures/generate_figure2_workflow.jl`：生成主算法流程图。
- `scripts/figures/generate_figure8_source_path.jl`：生成 100 节点 source 1 首次接收传播树。
- `scripts/figures/generate_large_scale_age_figures.jl`：生成大规模 AoI 统计图。

### 9.2 Fig. 3、Fig. 4、Fig. 5 与历史对比图

- `scripts/processing/CMP_sys_inf_nc_inf_poolA.jl`：小图 `sys/sys-inf/nc/nc-inf` AoI 对比图。
- `scripts/processing/INF_small_graphs_use_C_e.jl`：小规模协议干扰结果处理与绘图。
- `scripts/processing/INF_large_graphs_sys_sysinf_euc_big.jl`：欧氏图中规模 AoI 增益图。
- `scripts/processing/INF_large_graphs_sys_sysinf_euc_big_T.jl`：欧氏图中规模周期项增益图。
- `scripts/processing/INF_large_graphs_sys_sysinf_reny_big.jl`：Erdős-Rényi 图中规模 AoI 增益图。
- `scripts/processing/INF_large_graphs_sys_sysinf_reny_big_T.jl`：Erdős-Rényi 图中规模周期项增益图。
- `scripts/processing/small_graphs.jl`：原始小图结果处理。
- `scripts/processing/large_graphs.jl`：原始大图结果处理。
- `scripts/processing/inflight.jl`：in-flight 数量相关图表处理。
- `scripts/processing/plot_algorithm_comparison.py`：Python 版本算法对比绘图。

### 9.3 Fig. 6、Fig. 7、可扩展性和分布式图

- `scripts/processing/plot_algorithm_convergence_fig6.jl`：生成图 6 收敛曲线。
- `scripts/processing/finalize_scalability_statistics.jl`：聚合可扩展性统计并生成相关图表。
- `paper_anonymized_2026_7/reproducibility/plot_algorithm_convergence_fig6.jl`：论文 `_7` 版本图 6 的可复现实验绘图副本。
- `paper_anonymized_2026_7/reproducibility/finalize_scalability_statistics.jl`：论文 `_7` 版本可扩展性统计副本。

### 9.4 主要输出图文件

- `plots/five_node_pan_comparison.pdf`：五节点 pan 示例图。
- `plots/algorithm_workflow.pdf`：算法流程图。
- `plots/algorithm_convergence.pdf`：100 节点收敛图。
- `plots/algorithm_decomposition.pdf`：周期项和延迟项分解图。
- `plots/best_age_first_source_1.pdf`：source 1 首次接收传播树。
- `plots/scalability_statistics.pdf`：100/400/700/1000 可扩展性统计图。
- `plots/distributed_scalability.pdf`：分布式可扩展性图。
- `plots/*.eps`：部分 Fig. 3--Fig. 5 的 EPS 版本，供期刊矢量图提交使用。

## 10. 论文与可复现实验包

- `paper_anonymized_2026_7/myJCIN.tex`：当前主论文 LaTeX 源文件。
- `paper_anonymized_2026_7/myJCIN.pdf`：当前编译 PDF。
- `paper_anonymized_2026_7/reproducibility/`：论文 `_7` 版本的脚本、结果表和调度序列副本。
- `paper_anonymized_2026_7/figures/`：论文内部引用的部分图文件。
- `paper_anonymized_2026_7/模板包/`：目标期刊模板文件。
- `paper_anonymized_2026-6/`：上一版论文和扩展讨论，部分分布式材料从该版本迁移到 `_7`。
- `docs/ieee_uncoded_all_to_all_aoi.tex`：无编码 all-to-all AoI 建模说明文档。
- `docs/ieee_uncoded_all_to_all_aoi.pdf`：上述文档编译结果。

## 11. 文献与外部参考目录

- `sys文献/`：systematic scheduling、global interference、fundamental bounds 等文献。
- `混合线性规划+AoI/`：MILP、column generation、branch-and-price、Lagrangian decomposition 和 AoI 优化相关文献。
- `paper_anonymized_2026_7/局部贪婪调度/`：local greedy scheduling、GNN scheduling、分布式 shortest path 等分布式部署参考文献。
- `JCIN LaTex Template/` 与 `JCIN LaTex Template.zip`：目标期刊模板资料。

## 12. 数据目录约定

- `data/exp_raw/`：原始小图 `.jld2` 结果，数量较多，不在本文档逐一展开。
- `data/exp_raw_sys_inf_e/`：协议干扰增量 A* 或相关历史实验结果池。
- `data/exp_res/`：当前主要结果目录，存放 `.csv`、`.txt`、`.g6`、调度序列和验证输出。
- `paper_anonymized_2026_7/reproducibility/`：投稿论文使用的冻结副本，优先用于复现实验和审稿材料。
- `tmp/`、`tmp_inf/`、`tmp_inf_e/`：临时文件目录，可由脚本重新生成，不应作为论文结果的唯一来源。

## 13. 推荐复现顺序

1. 小图语义检查：运行 `scripts/evaluations/validate_pan5_inf_astar_v8.jl`。
2. 小图精确验证：运行 `scripts/evaluations/small_branch_price_cut_validation.jl`。
3. 100 节点集中式策略比较：运行 `scripts/evaluations/interference_algorithm_comparison.jl`。
4. 100 节点调度约束回放：运行 `scripts/evaluations/validate_saved_schedule_constraints.jl`。
5. 大规模统计：运行 `scripts/evaluations/scalability_statistics_ablation.jl`。
6. 大规模分布式部署：运行 `scripts/evaluations/distributed_scalability_experiment.jl`。
7. 论文图表：运行 `scripts/figures/*.jl` 和 `scripts/processing/*.jl` 中对应图的生成脚本。

对于耗时实验，建议优先使用 `data/exp_res` 和 `paper_anonymized_2026_7/reproducibility` 中已有结果；需要重新跑服务器批量实验时，再使用第 7 节中的 batch 入口。
