using Graphs

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



if length(m_nodes_metas) != length(m_vertices)
  println(length(m_vertices))
  println(length(m_nodes_metas))
  println("FFFFFF")
  
end
if length(nbd_nodes_metas) != length(nbd_vertices) 

  println("EEEE")
end



#### HERE have all the meta data and vertices

####### 