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


function apsp(adj_list :: Dict{Int, Set{Int}}, nodes :: Vector{Int})
  all_pairs_distances = Dict{Tuple{Int, Int}, Int}()  # Map (n1, n2) -> distance

  for source in nodes
      # Run BFS from the source node
      dist_from_source = bfs(source, adj_list)

      # Add distances to the dictionary
      for (target, distance) in dist_from_source
          all_pairs_distances[(source, target)] = distance
      end
  end

  return all_pairs_distances
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


# this will not include 0s and infs
function get_mean_path_length(apsp_list)
  total_pairs :: Int = 0
  total_dist :: Int = 0



  for ((n1, n2), dist) in apsp_list
    if dist != 0 && dist != 1000000000
      total_pairs += 1 
      total_dist += dist
    end 

  end

  average_dist :: Float64 = total_dist / total_pairs

  return average_dist


end 



function comp_global_cc(vertices, adj_dict) :: Float64
  c = 0
  pathlengths2 = 0
  for v in vertices
    neighs = adj_dict[v]
    for i in neighs, j in neighs
      i == j && continue
      if i in adj_dict[j]
        c += 1
      end
    end
    k = length(neighs)
    pathlengths2 += k * (k - 1)
  end
  return c / pathlengths2
end

function apsp_built_in(g)
  n = nv(g) 

  all_pairs_distances = Dict{Tuple{Int, Int}, Int}()

  for i in 1:n


    # Use Dijkstra's algorithm from vertex `i`
    sp = dijkstra_shortest_paths(g, i)  # Shortest paths from vertex `i`
    
    # Extract the distances and store them in the dictionary
    for j in 1:n
      # Only store valid distances (not Inf or unreachable)
      if i != j && sp.dists[j] != Inf
        all_pairs_distances[(i, j)] = sp.dists[j]
      end
    end
  end
  return all_pairs_distances
end 

function create_graph(edgeSet, length)
  g = SimpleGraph(length)

  for (u, v) in edgeSet
    add_edge!(g, u, v)
  end
  return g
end 





###### begin question #######


dataset = "social-online-small-2.txt"

edges :: Set{Tuple{Int, Int}} = get_edges(dataset)
adj_list :: Dict{Int, Set{Int}} = get_adj_list(edges)
vertices :: Vector{Int} = get_vertices(adj_list)

config_model_cc :: Vector{Float64} = []

config_model_mean :: Vector{Float64} = []

apsp_list = apsp_built_in(create_graph(edges, length(vertices)))
orig_mean_path_len = get_mean_path_length(apsp_list) 
orig_cc = comp_global_cc(vertices, adj_list)



# need 20m swaps for config model, need ~1000 random graphs

for i in 1:1000
  global edges 
  global adj_list
  global count
  global apsp_list
  if i % 10 == 0 
    println(i)
  end
  (edges, adj_list) = double_edge_swap(edges, adj_list, 20*length(edges))
  apsp_list = apsp_built_in(create_graph(edges, length(vertices)))
  mean_path_len = get_mean_path_length(apsp_list)
  cc = comp_global_cc(vertices, adj_list)
  push!(config_model_cc, cc)
  push!(config_model_mean, mean_path_len)
end 

using Plots

gr()
histogram(config_model_cc, bins=30, title="Friendship Network vs Null - C", xlabel="Clustering Coeficent (C)", ylabel="Frequency", legend=true, label="configuration model")
vline!([orig_cc], label = "network", linestyle = :solid, color = :green)
savefig(current(), "cc-online.pdf")

gr()
histogram(config_model_mean, bins=30, title="Friendship Network vs Null - ⟨l⟩", xlabel="Mean path length ⟨l⟩", ylabel="Frequency", legend=true, label="configuration model")
vline!([orig_mean_path_len], label = "network", linestyle = :solid, color = :green)
savefig(current(), "mp-online.pdf")

## https://icon.colorado.edu/networks 
## M. Fire, and R. Puzis, "Organization mining using online social networks." Networks and Spatial Economics 16(2), 545-578 (2016)
# small-2 = S1 *** this is the one I used




########## end ##############



###### begin other 

