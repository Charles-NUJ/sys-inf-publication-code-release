using DrWatson
@quickactivate "age-optimal-multisource-flooding"
using Graphs
include(joinpath(srcdir(), "solvers", "a_star_src_inf_v3.jl"))

graph = SimpleGraph(5)
for edge in [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)]
    add_edge!(graph, edge...)
end
connectivity = Matrix{Int8}(adjacency_matrix(graph))
solution_knowledge = Dict{String,Any}(
    "dissemination_delays" => zeros(5, 5),
    "t_star" => 12,
)

for epsilon in (1.0, 0.0)
    solution = a_star_inf_search(connectivity, 1, 2, 2,
                                 deepcopy(solution_knowledge), epsilon)
    first_decodings = copy(solution.decoding_levels)
    first_decodings[I(5)] .= 1
    delay_sum = sum(first_decodings[s, j] - solution.first_encodings[s]
                    for s in 1:5, j in 1:5 if s != j)
    period = length(solution.schedule)
    aoi = period / 2 + delay_sum / 20
    println("epsilon=$epsilon")
    println("schedule=$(solution.schedule)")
    println("payload_schedule=$(solution.payload_schedule)")
    println("T=$period delay_sum=$delay_sum AoI=$aoi")
end
