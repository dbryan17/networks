using Graphs
using Random
using Plots

dir_path = "./data/"

function get_edges(filename) :: Set{Tuple{Int, Int}}
  edge_set = Set()
  open(dir_path * filename) do f 
    for line in eachline(f)
      nodes = split(line)
      node1, node2 = parse(Int, nodes[1]), parse(Int, nodes[2])
      if node1 == 0 
        node1 = 1421
      end 
      if node2 == 0 
        node2 = 1421
      end
      edge = tuple(min(node1, node2), max(node1, node2))
      push!(edge_set, edge)
    end
  end
  return edge_set
end



# from edge set
function get_adj_list(edges :: Set{Tuple{Int, Int}}) :: Dict{Int, Set{Int}}
  adj_list = Dict{Int, Set{Int}}()
  for (n1, n2) in edges
    if (haskey(adj_list, n1))
      push!(adj_list[n1], n2)
    else
      adj_list[n1] = Set(n2)
    end

    if (haskey(adj_list, n2))
      push!(adj_list[n2], n1)
    else
      adj_list[n2] = Set(n1)
    end
  end
  return adj_list
end


# from adj dict
function get_vertices(adj_list :: Dict{Int, Set{Int}}) :: Vector{Int}
  keyss = keys(adj_list)
  vertices :: Vector{Int} = []
  for key in keyss
    push!(vertices, key)
  end
  return vertices
end


function remove_some_edges(edge_dict :: Dict{Int, Set{Int}}, edges :: Set{Tuple{Int, Int}}, frac_to_observe :: Float64) :: Tuple{Dict{Int, Set{Int}}, Set{Tuple{Int, Int}}, Set{Tuple{Int, Int}}}
  fration_to_remove = 1 - frac_to_observe
  num_to_remove :: Int = round(fration_to_remove * length(edges))
  edge_dict = deepcopy(edge_dict)
  edges = deepcopy(edges)
  rmed_edges = Set()

  for i in 1:num_to_remove
    rand_edge = rand(edges)
    pop!(edges, rand_edge)
    (n1, n2) = rand_edge
    pop!(edge_dict[n1], n2)
    pop!(edge_dict[n2], n1)
    push!(rmed_edges, rand_edge)
  end

  return (edge_dict, edges, rmed_edges)
end

function add_meta_data_nbd(filename) :: Dict{Int, Int}
  # returns a dict which is vertex -> metavar (1 or 2)
  node_meta :: Dict{Int, Int} = Dict()
  open(dir_path * filename) do f 
    for line in eachline(f)
      all_stuff = split(line, ",")
      node_meta[parse(Int,all_stuff[1])] = parse(Int, all_stuff[4])
    end
  end
  return node_meta
end

function add_meta_data_m(filename) :: Dict{Int, Int}
  node_meta :: Dict{Int, Int} = Dict()
  open(dir_path * filename) do f 
    for (i, line) in enumerate(eachline(f))
      node_meta[i] = parse(Int, line)
    end 
  end
  return node_meta
end

######## 
# nbd is 0-indexed, m is 1-indexed
nbd_edges_f = "nbd-edges.txt"
nbd_meta_f = "nbd-nodes.txt"

m_edges_f = "m-edges.txt"
m_meta_f = "m-meta.txt"

# 0-indexed
nbd_edges :: Set{Tuple{Int, Int}} = get_edges(nbd_edges_f)
nbd_adj_list :: Dict{Int, Set{Int}} = get_adj_list(nbd_edges)
nbd_vertices :: Vector{Int} = get_vertices(nbd_adj_list)

nbd_nodes_metas :: Dict{Int, Int} = add_meta_data_nbd(nbd_meta_f)



m_edges :: Set{Tuple{Int, Int}} = get_edges(m_edges_f)
m_adj_list :: Dict{Int, Set{Int}} = get_adj_list(m_edges)
m_vertices :: Vector{Int} = get_vertices(m_adj_list)

m_nodes_metas :: Dict{Int, Int} = add_meta_data_m(m_meta_f)

# add vertices with no edges
for (node, meta) in m_nodes_metas
  if !(node in m_vertices)
    push!(m_vertices, node)
    m_adj_list[node] = Set()
  end
end

# keep meta data just to add nodes with no edges



if length(m_nodes_metas) != length(m_vertices)
  println(length(m_vertices))
  println(length(m_nodes_metas))
  println("bad")
  
end
if length(nbd_nodes_metas) != length(nbd_vertices) 

  println("bad")
end


function jc(n1_neigh :: Set{Int}, n2_neigh :: Set{Int}):: Float64
  inter = intersect(n1_neigh, n2_neigh)
  un = union(n1_neigh, n2_neigh)

  return length(inter) / length(un)

end 

function dp(n1_neigh :: Set{Int}, n2_neigh :: Set{Int})
  return length(n1_neigh) * length(n2_neigh)
end

function create_graph(edgeSet, length)
  g = SimpleGraph(length)

  for (u, v) in edgeSet
    add_edge!(g, u, v)
  end
  return g
end 

function shortest_path_preder(length) 
  if length == 0
    # print("INFFF")
    return 0.0
  else
    return 1/length
  end
end


function get_edges_to_pred(edge_dict_obs :: Dict{Int, Set{Int}}, edges_obs :: Set{Tuple{Int, Int}}) :: Set{Tuple{Int, Int}}

  to_predict :: Set{Tuple{Int, Int}} = Set()
  for n1 in keys(edge_dict_obs)
    for n2 in keys(edge_dict_obs)
      if n1 != n2
        poss_edge = tuple(min(n1, n2), max(n1, n2))
        if !(poss_edge in edges_obs)
          push!(to_predict, poss_edge)
        end
      end
    end
  end
  return to_predict
end 





function make_predictions(edge_dict_obs :: Dict{Int, Set{Int}}, edges_obs :: Set{Tuple{Int, Int}}, edges_missing :: Set{Tuple{Int, Int}}) :: Vector{Vector{Float64}}
  # we make a prediction for each possible missing edge

  # this will be for a partially observed graph

  # edges_missing is the edges which are missing from this partially observed graph

  edges_to_predict :: Set{Tuple{Int, Int}} = get_edges_to_pred(edge_dict_obs, edges_obs)

  g = create_graph(edges_obs, length(edge_dict_obs))

  # keys are nodes, vertices are the weird thing that is returned form the built in shortest path alg

  shortest_path_from = Dict()
  for (i, j) in edges_to_predict
    if !(i in keys(shortest_path_from))
      shortest_path_from[i] = dijkstra_shortest_paths(g, i).pathcounts
    end
    if !(j in keys(shortest_path_from))
      shortest_path_from[j] = dijkstra_shortest_paths(g, j).pathcounts
    end
  end

  # i j tk JC DP SP TPR FPR 

  score_table :: Vector{Vector{Float64}}  = []
 
  for (i, j) in edges_to_predict

    # do predictions
    jc_score = jc(edge_dict_obs[i], edge_dict_obs[j]) * (rand() * 0.0001)
    dp_score = dp(edge_dict_obs[i], edge_dict_obs[j]) * (rand() * 0.0001)
    shortest_path = shortest_path_from[i][j]
    sp_score = shortest_path_preder(shortest_path) * (rand() * 0.0001)
    tau = 0
    if (i,j) in edges_missing
      tau = 1
    end
    push!(score_table, [i, j, tau, jc_score, dp_score, sp_score, -1, -1])
  end

  # TODO do real calculation of TPR and FPR

  return score_table

end



function two_a(edge_dict :: Dict{Int, Set{Int}}, edges :: Set{Tuple{Int, Int}})

  # function remove_some_edges(edge_dict :: Dict{Int, Set{Int}}, edges :: Set{Tuple{Int, Int}}, frac_to_observe :: Float64) :: Tuple{Dict{Int, Set{Int}}, Set{Tuple{Int, Int}}}

  # need full ROC curve for the three predictors for one of the 
  
  # TODO question for OH... do we need the do the thing where we break ties randomly and compute average? 
  # OR can we just use epsilon and call it good

  frac_to_obs = 0.05

  while frac_to_obs < 1
    (edge_dict_obs, edges_obs, rmed_edges) = remove_some_edges(edge_dict, edges, frac_to_obs)
    println(frac_to_obs)
    score_table = make_predictions(edge_dict_obs, edges_obs, rmed_edges)
    frac_to_obs += 0.05

  end
end

####### here start problem #########

two_a(nbd_adj_list, nbd_edges)
two_a(m_adj_list, m_edges)
####################################



# TODO add to discussion
# assuming if we did not observe any attrs, we make no predictions 
# the baseline doesn't work, and we could just pick randomly one of the attrs
# but I'm going on the assumption that we wouldn't even know what the attrs are

# create plots
# x1 = collect(keys(m_b_to_acc))
# y1 = collect(values(m_b_to_acc))

# x2 = collect(keys(nbd_b_to_acc))
# y2 = collect(values(nbd_b_to_acc))

# x1_sorted, y1_sorted = first.(sort(collect(m_b_to_acc), by=first)), last.(sort(collect(m_b_to_acc), by=first))
# x2_sorted, y2_sorted = first.(sort(collect(nbd_b_to_acc), by=first)), last.(sort(collect(nbd_b_to_acc), by=first))

# plot(x1_sorted, y1_sorted, label="Malaria var DBLa Cys-PoLV groups", lw=2, color = :red,
#     xlabel="β - amount of randomness",
#     ylabel="ACC - accurary of predictions (correct / total)",
#     title = "Local smoothing accuracy as a function of randomness", 
#     )
# plot!(x2_sorted, y2_sorted, label="Board of Directors Gender", lw=2, color = :green)

# hline!([0.5386345661562933], label = "BD baseline ACC", linestyle = :dash, color = :green)
# hline!([0.36977580663985826], label = "M baseline ACC", linestyle = :dash, color = :red)




# savefig(current(), "test-b.pdf")

####### 
