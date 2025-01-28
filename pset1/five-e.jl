
######### end header #########

using LightGraphs


# import 
dir_path = "./facebook100txt/"

# is 100 which is correct
# print(length(names))
# EDGES ARE DUPLICATED

get_names = () -> 
  # get names of ones that acutally contain data
  filter(file -> endswith(file, ".txt") && 
                          !endswith(file, "_attr.txt") && 
                          !occursin("readme", file),
                 readdir(dir_path))


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

function get_all_edges() :: Dict{String, Set{Tuple{Int, Int}}}
  ret = Dict()
  for name in get_names()
    ret[name[1:end-4]] = get_edges(name)
    
  end
  return ret
end


######### end header #########

names = get_names()

# TODO fix me for whole thing
names = [names[2]]
# println(names)


vertices :: Dict{String, Set{Int}} = Dict()

for name in names
  name_nodes = Set()
  open(dir_path * name) do f 
    for line in eachline(f) 
      nodes = split(line)
      node1, node2 = parse(Int, nodes[1]), parse(Int, nodes[2])
      push!(name_nodes, node1)
      push!(name_nodes, node2)
    end
  end 
  short_name = name[1:end-4]
  vertices[short_name] = name_nodes
end

# TODO change to get all edges
# edges :: Dict{String, Set{Tuple{Int, Int}}} = get_all_edges()
edges :: Dict{String, Set{Tuple{Int, Int}}} = Dict()
edges["Amherst41"] = get_edges("Amherst41.txt")


# Function to perform DFS and find a connected component
function dfs(v, adj_list, visited, component)
  push!(component, v)
  visited[v] = true
  for neighbor in adj_list[v]
      if !get(visited, neighbor, false)
          dfs(neighbor, adj_list, visited, component)
      end
  end
end

function bfs(source::Int, adj_list::Dict{Int, Set{Int}})
  dist = Dict{Int, Int}()  # Dictionary to handle arbitrary node IDs

  # Initialize distances to Inf for all nodes
  for node in keys(adj_list)
      dist[node] = 10000000
  end
  dist[source] = 0.0  # Distance to itself is 0

  # BFS queue
  queue = [source]
  while !isempty(queue)
      current = popfirst!(queue)
      for neighbor in adj_list[current]
          if dist[neighbor] == 10000000  # Check if the neighbor hasn't been visited
              dist[neighbor] = dist[current] + 1
              push!(queue, neighbor)
          end
      end
  end
  return dist
end

function apsp(adj_list, nodes)
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

function dijkstra(source, adj_list, nodes)
  dist = Dict{Int, Float64}()
  prev = Dict{Int, Int}()
  unvisited = Set{Int}(nodes)  # Set of unvisited nodes

  # Initialize distances with infinity and previous nodes with none
  for v in nodes
      dist[v] = Inf  # Use Inf (as Float64)
      prev[v] = -1
  end
  dist[source] = 0.0  # Distance from the source to itself is 0

  while !isempty(unvisited)
      # println("Unvisited nodes remaining: ", length(unvisited))

      # Get the node in `unvisited` with the smallest distance
      u = argmin(v -> dist[v], unvisited)

      delete!(unvisited, u)  # Correctly remove node u from unvisited

      # Update the distances of the neighbors of u
      for v in adj_list[u]
          alt = dist[u] + 1  # Assuming edge weights are all 1 for simplicity
          if alt < dist[v]
              dist[v] = alt
              prev[v] = u
          end
      end
  end

  return dist
end


figure_max_total_size :: Dict{String, Tuple{Int, Int}} = Dict()
figure_avg_comp_size :: Dict{String, Tuple{Float64, Int}} = Dict()


for (name, edges) in edges

  println(name)

  ad_list :: Dict{Int, Set{Int}} = Dict()

  # init with nodes 
  for node in vertices[name]
    ad_list[node] = Set()

  end

  for (n1, n2) in edges
    push!(ad_list[n1], n2)
    push!(ad_list[n2], n1)
  end

  # Find all connected components
  visited = Dict{Int, Bool}()
  components = []

  for v in vertices[name]
    if !get(visited, v, false)
        component = []
        dfs(v, ad_list, visited, component)
        push!(components, component)
    end
  end

  max = 0
  sum = 0
  largest_comp = []
  for comp in components
    sum += length(comp)
    if max < length(comp)
      largest_comp = comp
      max = length(comp)
    end
  end

  new_adj_list :: Dict{Int, Set{Int}} = Dict()

  for n in largest_comp 
    new_adj_list[n] = ad_list[n]
  end 
  # this is working up to here


  count = length(largest_comp)

  all_pairs_distances = apsp(new_adj_list, largest_comp)

  # print(dist_matrix)

  # Find All-Pairs Shortest Path in the largest component
  # all_pairs_distances = Dict{Tuple{Int, Int}, Int}()
  # for u in largest_comp
  #     flag = false 
  #     if count % 100 == 0 
  #       println("count ")
  #       print(count)
  #       flag = true 
  #     end
  #     count -= 1
  #     dist = dijkstra(u, new_adj_list, largest_comp)
  #     for v in largest_comp
  #         # if flag
  #         #   println(dist[v])
  #         # end
  #         all_pairs_distances[(u, v)] = dist[v]
  #     end
  # end

  total_dist :: Int = 0

  max_dist :: Int = 0


  total_pairs :: Int = 0

  size_of_largest_comp :: Int = length(largest_comp)

  for ((n1, n2), dist) in all_pairs_distances
    if(n1 != n2)

      total_dist += dist
      total_pairs += 1

      if dist > max_dist
        max_dist = dist
      end
    end
  end

  average_dist :: Float64 = total_dist / total_pairs


  println("max_dist")
  println(max_dist)

  println("average_dist")
  println(average_dist)

  println("size_of_largest_comp")
  println(size_of_largest_comp)

  figure_avg_comp_size[name] = (average_dist, size_of_largest_comp)
  figure_max_total_size[name] = (max_dist, length(vertices[name]))

end



# figure_max_total_size :: Dict{String, Tuple{Int, Int}} = Dict()
# figure_avg_comp_size :: Dict{String, Tuple{Float64, Int}} = Dict()

max_fig_x = []
max_fig_y = []
for (name, (max, total)) in figure_max_total_size
  push!(max_fig_x, total)
  push!(max_fig_y, max)
end


avg_fig_x = []
avg_fig_y = []

for (name, (avg, largest)) in figure_avg_comp_size
  push!(avg_fig_x, largest)
  push!(avg_fig_y, avg)
end

using Plots
gr()
max_fig = scatter(max_fig_x, max_fig_y, title="Diameter as a function of total network size", xlabel="network size", ylabel="diameter", color=:blue, markersize=5)
avg_fig = scatter(avg_fig_x, avg_fig_y, title="Mean geodesic distance (of largest componenet) as a function of largest component size", xlabel="largest component size", ylabel="mean geodesic distance", color=:blue, markersize=5)
savefig(max_fig, "e-maxfig.pdf")
savefig(avg_fig, "e-avgfig.pdf")




