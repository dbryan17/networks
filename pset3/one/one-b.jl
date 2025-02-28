# Lots of copying from one-b.jl

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

function double_edge_swap(edges :: Set{Tuple{Int, Int}}, adj_dict :: Dict{Int, Set{Int}}, desired :: Int) :: Tuple{Set{Tuple{Int, Int}}, Dict{Int, Set{Int}}}
  # perform [desired] double edge swaps of the given simple graph with edge set and adj_list
  
  adj_dict = deepcopy(adj_dict)
  edges = deepcopy(edges)


  edge_list = collect(edges)
  edges_len = length(edge_list)
  successfull = 0




  while successfull < desired
    
    # get random edges
    idx1, idx2 = rand(1:edges_len), rand(1:edges_len)
    while idx1 == idx2
      idx2 = rand(1:edges_len)
    end

    (u, v) = edge_list[idx1]
    (x, y) = edge_list[idx2]

    # coin flip for diangnol or vertical
    vertSwap = rand(Bool)

    # could make it so if vert swap fails, try diagnol, but I think that makes this not random

    if vertSwap 
      if (u != x && v != y) && !(u in adj_dict[x]) && !(v in adj_dict[y]) # prevent self loops and multigraph
        successfull += 1

        # delete old from edge_list
        filter!(edge -> edge != (u,v), edge_list)
        filter!(edge -> edge != (x,y), edge_list)

        # delete old from dict
        delete!(adj_dict[u], v)
        delete!(adj_dict[v], u)
        delete!(adj_dict[x], y)
        delete!(adj_dict[y], x)

        # add to edge list
        push!(edge_list, (u,x))
        push!(edge_list, (v,y))

        # add to dict 
        push!(adj_dict[u], x)
        push!(adj_dict[x], u)
        push!(adj_dict[v], y)
        push!(adj_dict[y], v)
      end
    else 

      # check if diagnol works 
      if (u != y && v != x) && !(u in adj_dict[y]) && !(v in adj_dict[x]) # check
        successfull += 1

        # delete old from edge_list
        filter!(edge -> edge != (u,v), edge_list)
        filter!(edge -> edge != (x,y), edge_list)

        # delete old from dict
        delete!(adj_dict[u], v)
        delete!(adj_dict[v], u)
        delete!(adj_dict[x], y)
        delete!(adj_dict[y], x)

        # add to edge list
        push!(edge_list, (u,y))
        push!(edge_list, (v,x))

        # add to dict 
        push!(adj_dict[u], y)
        push!(adj_dict[y], u)
        push!(adj_dict[v], x)
        push!(adj_dict[x], v)

      end
    end
  end


  # remake edges set
  edges = Set(edge_list)

  return (edges, adj_dict)

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

  curr_a = .01

  num_iters = 100

  while curr_a < 1.0
    total_avg = 0
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


function calculate_b(r :: Int, m :: Int) :: Float64
  return r / (2*m)
end


function oneb(nodes_meta_truth :: Dict{Int, Int}, edge_dict :: Dict{Int, Set{Int}}, edges :: Set{Tuple{Int, Int}}) :: Dict{Float64, Float64}

  num_swaps :: Int = 0  

  edges_len = length(edges)

  num_iters = 1000

  num_bs = 100

  max_num_swaps = 2 * edges_len

  # just evenly spaced
  swap_intervals = round.(Int, range(10, max_num_swaps; length=num_bs))

  b_to_acc :: Dict{Float64, Float64} = Dict()


  total_acc = 0
  for i in 1:num_iters
    nodes_metas_obs = remove_some_metas(nodes_meta_truth, .8)
    acc = local_smooth(nodes_metas_obs, nodes_meta_truth, edge_dict)
    total_acc += acc
  end
  average_acc = total_acc / num_iters

  # uncomment for non log 10
  b_to_acc[calculate_b(0, edges_len)] = average_acc


  for num_swaps in swap_intervals

    println("curr num swaps out of 7710 if diff is 78")
    println("out of 5368 if diff is 54")
    println(num_swaps)

    total_acc = 0
    for i in 1:num_iters
      (rand_edges, rand_adj_list) = double_edge_swap(deepcopy(edges), deepcopy(edge_dict), num_swaps)
      nodes_metas_obs = remove_some_metas(nodes_meta_truth, .8)
      acc = local_smooth(nodes_metas_obs, nodes_meta_truth, rand_adj_list)
      total_acc += acc
    end

    average_acc = total_acc / num_iters

    b_to_acc[calculate_b(num_swaps, edges_len)] = average_acc

    num_swaps += round(max_num_swaps / num_bs)

  end

  return b_to_acc

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


# β = r/2m   r is number of double swaps tha


nbd_b_to_acc = oneb(nbd_nodes_metas, nbd_adj_list, nbd_edges)

m_b_to_acc :: Dict{Float64, Float64} = oneb(m_nodes_metas, m_adj_list, m_edges)
# 0.36977580663985826 --- m 
# 0.5386345661562933 --- nbd


# create plots
x1 = collect(keys(m_b_to_acc))
y1 = collect(values(m_b_to_acc))

x2 = collect(keys(nbd_b_to_acc))
y2 = collect(values(nbd_b_to_acc))

x1_sorted, y1_sorted = first.(sort(collect(m_b_to_acc), by=first)), last.(sort(collect(m_b_to_acc), by=first))
x2_sorted, y2_sorted = first.(sort(collect(nbd_b_to_acc), by=first)), last.(sort(collect(nbd_b_to_acc), by=first))

plot(x1_sorted, y1_sorted, label="Malaria var DBLa Cys-PoLV groups", lw=2, color = :red,
    xlabel="β - amount of randomness",
    ylabel="ACC - accurary of predictions (correct / total)",
    title = "Local smoothing accuracy as a function of randomness", 
    )
plot!(x2_sorted, y2_sorted, label="Board of Directors Gender", lw=2, color = :green)

hline!([0.5386345661562933], label = "BD baseline ACC", linestyle = :dash, color = :green)
hline!([0.36977580663985826], label = "M baseline ACC", linestyle = :dash, color = :red)




savefig(current(), "one-b.pdf")

####### 
