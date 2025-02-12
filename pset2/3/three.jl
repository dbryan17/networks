# question 3 
# r = the number of double edge swaps we have applied to some input graph G 
# carry out numverical experiment to answer the following question: 
# as a function of r, how does the clustering coeffiecent C, and mean path length <l> relax 
# onto those of the corresponding configuration model
# as references, overlay in your plots horizontal lines for the C and <l> when you have r = 20m swaps
# comment on how "random" Berkley's social network was to begin with, in in what ways, 
# and on the rate at which randomization destroys the empricial patterns

# import 
dir_path = "../data/"


# from file names
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
  keys = keys(adj_list)
  vertices :: Vector{Int} = []
  for key in keys
    push!(vertices, key)
  end
  return vertices
end



function dfs_iterative(v, adj_list, visited, component)
  stack = [v]  # Use an explicit stack
  while !isempty(stack)
      node = pop!(stack)
      if !get(visited, node, false)
          push!(component, node)
          visited[node] = true
          for neighbor in adj_list[node]
              if !get(visited, neighbor, false)
                  push!(stack, neighbor)  # Push unvisited neighbors
              end
          end
      end
  end
end

function double_edge_swap() 
  # given a simple graph, perform one double edge swap 


end



###### begin question #######

dataset = "Amherst41.txt"


# graphs will be adjacney lists

edges :: Set{Tuple{Int, Int}} = get_edges(dataset)
adj_list :: Dict{Int, Set{Int}} = get_adj_list(edges)
vertices :: Vector{Int} = get_vertices(adj_list)



########## end ##############


