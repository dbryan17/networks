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


function remove_some_metas(nodes_metas :: Dict{Int, Int}, fraction_obs :: Float64) :: Dict{Int, Int}

  # fraction obs is the fraction we want to have observered

  fraction_to_remove = 1 - fraction_obs
  num_to_remove :: Int = round(fraction_to_remove * length(nodes_metas))

  old_dict = deepcopy(nodes_metas)

  rand_keys = shuffle(collect(keys(old_dict)))

  keys_to_rm = rand_keys[1: num_to_remove]

  for k in keys_to_rm
    pop!(old_dict, k)
  end

  return old_dict
end


function predict_meta(node :: Int, edge_dict :: Dict{Int, Set{Int}}, nodes_metas_obs :: Dict{Int, Int}, all_attrs :: Vector{Int}) :: Int 


  # return rand(all_attrs)
  edges = edge_dict[node]

  

  counts :: Vector{Int} = []
  for edge in edges
    if haskey(nodes_metas_obs, edge)
      val = nodes_metas_obs[edge]
      push!(counts, val)
    end
  end

  if length(counts) == 0
    # pull from distrubtion of all_attrs... which is the same as a random draw from the array
    return rand(all_attrs)
  end

  counts_dict = Dict()

  for meta in counts
    if haskey(counts_dict, meta)
      counts_dict[meta] = counts_dict[meta] + 1
    else 
      counts_dict[meta] = 1
    end
  end

  maxes :: Vector{Int} = []
  max = -1
  

  for (meta_var, count) in counts_dict
    if count > max
      max = count
    end
  end

  for (meta_var, count) in counts_dict
    if max == count
      push!(maxes, meta_var)
    end
  end

  # now return a random meta variable from all the maxes
  return rand(maxes)

end


function local_smooth(nodes_metas_obs :: Dict{Int, Int}, nodes_metas_truth :: Dict{Int, Int}, edges_dict :: Dict{Int, Set{Int}}) :: Float64
  # get all the attrs for baseline
  all_attrs :: Vector{Int} = []
  for (_, val) in nodes_metas_obs
    push!(all_attrs, val)
  end

  correct_pred = 0
  number_pred = 0

  for (node, _) in nodes_metas_truth

    if !(haskey(nodes_metas_obs, node))
      predicted_meta_var = predict_meta(node, edges_dict, nodes_metas_obs, all_attrs)
      real_meta_var = nodes_metas_truth[node]
      if predicted_meta_var == real_meta_var
        correct_pred += 1
      end
      number_pred += 1
    end
  end
  return correct_pred / number_pred
end


function smooth_all(nodes_metas_truth :: Dict{Int, Int}, edges_dict :: Dict{Int, Set{Int}}) :: Dict{Float64, Float64}


  a_to_acc = Dict()

  # if we observed none, we will assume 0 percent accuracy
  # a_to_acc[0.0] = 0.0
  # we will assume 100 perfent accurary if we observe all of them, even tho there is
  # no inputs so this is really undef
  # a_to_acc[1.0] = 1.0 

  curr_a = .01

  num_iters = 1000

  while curr_a < 1.0
    total_avg = 0
    println("curr a ")
    println(curr_a)
    for i in 1:num_iters
      nodes_metas_obs = remove_some_metas(nodes_metas_truth, curr_a)
      total_avg += local_smooth(nodes_metas_obs, nodes_metas_truth, edges_dict)
    end
    full_avg = total_avg / num_iters
    a_to_acc[curr_a] = full_avg

    curr_a += .01
  end

  return a_to_acc

end


function compute_baseline_acc(nodes_metas_truth :: Dict{Int, Int}) :: Float64

  attrs_to_total :: Dict{Int, Int} = Dict()

  n = length(nodes_metas_truth)

  for attr in values(nodes_metas_truth)
    attrs_to_total[attr] = 0
  end

  for (node, attr) in nodes_metas_truth
    attrs_to_total[attr] = attrs_to_total[attr] + 1
  end

  sum = 0.0
  for (attr, total) in attrs_to_total
    sum += (total / n)^2
  end


  return sum

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


println(compute_baseline_acc(nbd_nodes_metas))
println(compute_baseline_acc(m_nodes_metas))


#### HERE have all the meta data and vertices

m_a_to_acc :: Dict{Float64, Float64} = smooth_all(m_nodes_metas, m_adj_list)
nbd_a_to_acc = smooth_all(nbd_nodes_metas, nbd_adj_list)




# TODO add to discussion
# assuming if we did not observe any attrs, we make no predictions 
# the baseline doesn't work, and we could just pick randomly one of the attrs
# but I'm going on the assumption that we wouldn't even know what the attrs are

# create plots
x1 = collect(keys(m_a_to_acc))
y1 = collect(values(m_a_to_acc))

x2 = collect(keys(nbd_a_to_acc))
y2 = collect(values(nbd_a_to_acc))

x1_sorted, y1_sorted = first.(sort(collect(m_a_to_acc), by=first)), last.(sort(collect(m_a_to_acc), by=first))
x2_sorted, y2_sorted = first.(sort(collect(nbd_a_to_acc), by=first)), last.(sort(collect(nbd_a_to_acc), by=first))

plot(x1_sorted, y1_sorted, label="Malaria var DBLa Cys-PoLV groups", lw=2, 
    xlabel="α - fraction of observed node labels compared to total",
    ylabel="ACC - accurary of predictions (correct / total)",
    title = "Effectiveness of local smoothing"
    )
plot!(x2_sorted, y2_sorted, label="Board of Directors Gender", lw=2)


savefig(current(), "one-a.pdf")

####### 
