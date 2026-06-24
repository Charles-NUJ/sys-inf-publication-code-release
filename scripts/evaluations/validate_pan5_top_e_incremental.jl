using DrWatson
quickactivate("age-optimal-multisource-flooding")

include(joinpath(scriptsdir(), "evaluations", "top_e_src_inf_graph_analysis_inf.jl"))

config = Dict{String,Any}(
    "piggyback_enable" => false,
    "max_coding_degree" => 1,
    "graph_id" => 1,
    "max_inflight" => 2,
    "cutoff" => "",
    "eps" => 0.0,
    "graph_size" => 5,
    "heuristic" => true,
    "graph_type" => "pan",
    "inf" => true,
    "max_inference_simu" => 2,
    "regenerate_random" => false,
)

run_experiment(5, config)
