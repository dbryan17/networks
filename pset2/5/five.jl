

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



#########################
dataset = "medici_network.gml"

id_to_name :: Dict{Int, String} = get_mapping(dataset)
edges :: Set{Tuple{Int, Int}} = get_edges(dataset)
adj_list :: Dict{Int, Set{Int}} = get_adj_list(edges)
adj_list[11] = Set()
vertices :: Vector{Int} = get_vertices(adj_list)


hcs :: Dict{Int, Float64} = Dict()

for node in vertices
  hcs[node] = comp_one_harmoic_cent(adj_list, node)
end

using Plots, StatsPlots
gr()

plot(xlabel = "Family Index", ylabel = "Harmonic Centrality",
     title = "Florentine Families Network Harmonic Centrality", legend = true, grid = true
    )


hc_other_i :: Vector{Int} = []
hc_other_hc :: Vector{Float64} = []

for (i, hc) in hcs
  if i == 8
    scatter!([i], [hc], label = "Medici family", color = :green, markersize = 5)

  else 
    push!(hc_other_i, i)
    push!(hc_other_hc, hc)

  end
end



scatter!(hc_other_i, hc_other_hc, label = "other family", color = :blue, markersize = 5)


savefig(current(), "hc.pdf")







#########################