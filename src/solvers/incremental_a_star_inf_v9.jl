#
# Parameterized interface for the AoI-aligned incremental interference solver.
#
# v9 keeps the v8 repair/search logic and exposes the two K-type controls used
# in experiments:
#   max_inflight       -- maximum number of unfinished source packets.
#   max_inference_simu -- maximum number of concurrent protocol-feasible
#                         broadcasters in one slot.
#

using Distributed
include(joinpath(@__DIR__, "incremental_a_star_inf_v8.jl"))

@everywhere begin
	using DrWatson
	quickactivate("age-optimal-multisource-flooding")

	function _v9_solution_knowledge(connectivity, solution_knowledge)
		if solution_knowledge === nothing
			node_count = size(connectivity, 1)
			return Dict{String, Any}(
				"dissemination_delays" => zeros(node_count, node_count),
				"t_star" => 0,
			)
		end
		return solution_knowledge
	end

	function incremental_a_star_inf_v9_search(connectivity;
											 max_coding_degree = 1,
											 max_inflight = 2,
											 max_inference_simu = 2,
											 solution_knowledge = nothing,
											 star_epsilon = 0.0,
											 repair_starts = 500,
											 repair_seed = 10_000,
											 repair_noise_scale = 0.65,
											 maximum_slots = 1_000_000,
											 run_base = true)
		@info "e-inf v9" max_inflight max_inference_simu
		return incremental_a_star_inf_v8_search(
			connectivity,
			max_coding_degree,
			max_inflight,
			max_inference_simu,
			_v9_solution_knowledge(connectivity, solution_knowledge),
			star_epsilon;
			repair_starts = repair_starts,
			repair_seed = repair_seed,
			repair_noise_scale = repair_noise_scale,
			maximum_slots = maximum_slots,
			run_base = run_base,
		)
	end

	function incremental_a_star_inf_v9_search(connectivity, max_coding_degree,
											 max_inflight, max_inference_simu,
											 solution_knowledge,
											 star_epsilon = 0.0;
											 repair_starts = 500,
											 repair_seed = 10_000,
											 repair_noise_scale = 0.65,
											 maximum_slots = 1_000_000,
											 run_base = true)
		return incremental_a_star_inf_v9_search(
			connectivity;
			max_coding_degree = max_coding_degree,
			max_inflight = max_inflight,
			max_inference_simu = max_inference_simu,
			solution_knowledge = solution_knowledge,
			star_epsilon = star_epsilon,
			repair_starts = repair_starts,
			repair_seed = repair_seed,
			repair_noise_scale = repair_noise_scale,
			maximum_slots = maximum_slots,
			run_base = run_base,
		)
	end
end
