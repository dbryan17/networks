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


function remove_some_edges(edge_dict :: Dict{Int, Set{Int}}, edges :: Set{Tuple{Int, Int}}, frac_to_observe :: Float64) :: Tuple{Dict{Int, Set{Int}}, Set{Tuple{Int, Int}}}
  fration_to_remove = 1 - frac_to_observe
  num_to_remove :: Int = round(fration_to_remove * length(edges))
  edge_dict = deepcopy(edge_dict)
  edges = deepcopy(edges)

  for i in 1:num_to_remove
    rand_edge = rand(edges)
    pop!(edges, rand_edge)
    (n1, n2) = rand_edge
    pop!(edge_dict[n1], n2)
    pop!(edge_dict[n2], n1)
  end

  return (edge_dict, edges)
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
  println("FFFFFF")
  
end
if length(nbd_nodes_metas) != length(nbd_vertices) 

  println("EEEE")
end



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
