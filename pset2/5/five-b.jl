

dir_path = "./data/"


function get_edges(filename) :: Set{Tuple{Int, Int}}
  edge_set = Set()
  currSrc = - 1
  open(dir_path * filename) do f 
    for line in eachline(f)
      if occursin("source", line)
        currSrc = parse(Int, split(line)[2])
      end
      if occursin("target", line)
        push!(edge_set, (currSrc, parse(Int, split(line)[2])))
      end
      # println(line)
      # nodes = split(line)

      # node1, node2 = parse(Int, nodes[1]), parse(Int, nodes[2])
      # edge = tuple(min(node1, node2), max(node1, node2))
      # push!(edge_set, edge)
    end
  end
  return edge_set
end


function get_mapping(filename) :: Dict{Int, String}
  currId = -1
  ret = Dict()
  open(dir_path * filename) do f
    for line in eachline(f)
      if occursin("id", line) && !occursin("Rid", line)
        currId = parse(Int, split(line)[2])
      end
      if occursin("label", line)
        ret[currId] = (split(line)[2])[2:end-1]
      end
    end
  end
  return ret
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


function bfs(source::Int, adj_list::Dict{Int, Set{Int}})
  dist = Dict{Int, Int}()  # Dictionary to handle arbitrary node IDs

  # Initialize distances to Inf for all nodes
  for node in keys(adj_list)
      dist[node] = 1000000000
  end
  dist[source] = 0  # Distance to itself is 0

  # BFS queue
  queue = [source]
  while !isempty(queue)
      current = popfirst!(queue)
      for neighbor in adj_list[current]
          if dist[neighbor] == 1000000000  # Check if the neighbor hasn't been visited
              dist[neighbor] = dist[current] + 1
              push!(queue, neighbor)
          end
      end
  end
  return dist
end

function comp_one_harmoic_cent(adj_dict :: Dict{Int, Set{Int}}, vertex :: Int) :: Float64
  lengths =  values(bfs(vertex, adj_list))

  sum :: Float64 = 0

  for length in lengths
    # non self and if it is inf, we will just add 0 so no effect
    if length != 0 && length != 1000000000
      sum += (1 / length)
    end
  end
  
  ci :: Float64 = (1 / (length(keys(adj_dict)) - 1)) * sum
  return ci 
end

function double_edge_swap(edges :: Set{Tuple{Int, Int}}, adj_dict :: Dict{Int, Set{Int}}, desired :: Int) :: Tuple{Set{Tuple{Int, Int}}, Dict{Int, Set{Int}}}
  # perform [desired] double edge swaps of the given simple graph with edge set and adj_list
  
  adj_dict = deepcopy(adj_dict)
  edges = deepcopy(edges)


  edge_list = collect(edges)
  edges_len = length(edge_list)
  successfull = 0




  while successfull < desired

    # println(length(edge_list))
    
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
        # println(length(edge_list))
        # println("1")
        filter!(edge -> edge != (u,v), edge_list)
        filter!(edge -> edge != (x,y), edge_list)
        # println(length(edge_list))
        # println("2")


        # delete old from dict
        delete!(adj_dict[u], v)
        delete!(adj_dict[v], u)
        delete!(adj_dict[x], y)
        delete!(adj_dict[y], x)

        # add to edge list
        push!(edge_list, (u,x))
        push!(edge_list, (v,y))
        # println(length(edge_list))
        # println("3")


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




#########################
dataset = "medici_network.gml"

id_to_name :: Dict{Int, String} = get_mapping(dataset)
edges :: Set{Tuple{Int, Int}} = get_edges(dataset)
adj_list :: Dict{Int, Set{Int}} = get_adj_list(edges)
adj_list[11] = Set()
vertices :: Vector{Int} = get_vertices(adj_list)


hcs :: Dict{Int, Float64} = Dict()

null_hcs :: Dict{Int, Vector{Float64}} = Dict()

null_hc_minus :: Dict{Int, Vector{Float64}} = Dict()

for node in vertices
  hcs[node] = comp_one_harmoic_cent(adj_list, node)
end

# init dict
for node in vertices
  null_hc_minus[node] = []
  null_hcs[node] = []
end

for i in 1:1000
  global edges
  global adj_list
  global count
  (edges, adj_list) = double_edge_swap(edges, adj_list, 20*length(edges))
  for node in vertices
    null_hc = comp_one_harmoic_cent(adj_list, node)
    push!(null_hcs[node], null_hc)
    push!(null_hc_minus[node], null_hc - hcs[node])
  end

end

using Plots, StatsPlots, CategoricalArrays
gr()
# ordered_nodes = []
ordered_null_hcs_minus :: Vector{Vector{Float64}} = []

ordered_null :: Vector{Vector{Float64}} = []
m_minus = -1
i = 1
for (node, hc_minus) in null_hcs
  global i 
  global m_minus
  if (node == 8)
    m_minus = i
  end
  i += 1
  
  # push!(ordered_nodes, node + 1)
  push!(ordered_null_hcs_minus, hc_minus)
end


ordered_nodes = [1:16]
# ordered_nodes = categorical(ordered_nodes)

println(ordered_nodes)

# println( length(ordered_nodes) == length(order_null_hc_minus) )

# println(ordered_null_hcs_minus)

# println(length(ordered_null_hcs_minus[2]))

# # Create the box plot
p = boxplot(ordered_nodes, ordered_null_hcs_minus, xlabel="Family Index", ylabel="Null HC - Actual", title="Florentine Families HC Null Minus", legend = false)
# colors = repeat(["#FF6347", "#4682B4", "#32CD32", "#FFD700", "#8A2BE2"], outer=length(ordered_nodes))


# # Add a vertical line at a specific x value (e.g., at the x-position of "C")
vline!(p, [m_minus], label="Medici", color=:green, linestyle=:dash)


savefig(p, "hc-minus-null.pdf")


# Example data: 16 x values with corresponding y-values
# x_values = [1:16]  # 16 categories
# y_values = [
#     [1.2, 1.4, 1.5, 2.1, 1.9],   # Values for "A"
#     [2.1, 2.4, 2.5, 2.2, 2.3],   # Values for "B"
#     [3.0, 3.2, 3.5, 3.8, 3.9],   # Values for "C"
#     [4.1, 4.3, 4.7, 4.9, 4.6],   # Values for "D"
#     [5.0, 5.2, 5.3, 5.5, 5.6],   # Values for "E"
#     [2.0, 2.1, 2.3, 2.5, 2.4],   # Values for "F"
#     [1.0, 1.2, 1.3, 1.5, 1.7],   # Values for "G"
#     [3.2, 3.4, 3.5, 3.7, 3.6],   # Values for "H"
#     [4.0, 4.2, 4.4, 4.6, 4.8],   # Values for "I"
#     [5.0, 5.3, 5.4, 5.6, 5.8],   # Values for "J"
#     [2.5, 2.6, 2.7, 2.9, 3.0],   # Values for "K"
#     [1.8, 2.0, 2.1, 2.3, 2.4],   # Values for "L"
#     [3.1, 3.3, 3.4, 3.5, 3.7],   # Values for "M"
#     [4.0, 4.1, 4.2, 4.3, 4.4],   # Values for "N"
#     [5.2, 5.3, 5.4, 5.5, 5.7],   # Values for "O"
#     [6.1, 6.3, 6.4, 6.5, 6.6]    # Values for "P"
# ]

# # Create the box plot
# test = boxplot(x_values, y_values, xlabel="Category", ylabel="Values", title="Box Plot with 16 Categories", legend=false)

# savefig(test, "test.pdf")





#########################