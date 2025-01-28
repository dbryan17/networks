# TODO figure out modules resuable stuff, etc... in julia


######### end header #########

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

### b ###


names = get_names()

# TODO fix me for whole thing
# names = [nameees[1], nameees[2]]
# println(names)

# get vertices 
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

# all edges 
# TODO fix me for whole thing
edges :: Dict{String, Set{Tuple{Int, Int}}} = get_all_edges()

# edges :: Dict{String, Set{Tuple{Int, Int}}} = Dict()
# edges["American75"] = edgeeees["American75"]
# edges["Amherst41"] = edgeeees["Amherst41"]
# println("edges")
# println(keys(edges))
# println("vertices")
# println(keys(vertices))
# println(edges["American75"])

# println("edges...")
# println(edges)

# println("vertices...")
# println(vertices)

# tuple is (vertex#, degree^2)
squared_degrees :: Dict{String, Set{Tuple{Int, Int}}} = Dict()

for (school, nodes) in vertices
  # school level
  println(school)
  print(haskey(edges, school))
  edge_set :: Set{Tuple{Int, Int}} = edges[school]
  degree_set :: Set{Tuple{Int, Int}} = Set()
  # key = node# value = degree
  degree_dict :: Dict{Int, Int} = Dict()
  for (n1, n2) in edge_set
    # edge level
    if haskey(degree_dict, n1)
      degree_dict[n1] += 1
    else 
      degree_dict[n1] = 1
    end
    if haskey(degree_dict, n2)
      degree_dict[n2] += 1
    else 
      degree_dict[n2] = 1
    end
  end
  # school level 
  for node in vertices[school]
    push!(degree_set, (node, (degree_dict[node])^2))
  end
  squared_degrees[school] = degree_set
end

edges_counts :: Dict{String, Int} = Dict()
for (k,v) in edges
  edges_counts[k] = length(v)
end

node_counts :: Dict{String, Int} = Dict()
for (k,v) in vertices
  node_counts[k] = length(v)
end

# <k_u>
average_degrees :: Dict{String, Float64} = Dict()
# school, 2m/n
for (school, node_count) in node_counts
  average_degrees[school] = 2*edges_counts[school] / node_count
end

average_squared_degrees :: Dict{String, Float64} = Dict()
for (school, squared_degree) in squared_degrees
  sum_squared_deg = 0
  for (node, sq_deg) in squared_degree
    sum_squared_deg += sq_deg
  end
  average_squared_degrees[school] = sum_squared_deg / node_counts[school]
end

# <k_v>
average_adjacent_degrees :: Dict{String, Float64} = Dict()
for (school, average_squared_degree) in average_squared_degrees
  average_adjacent_degrees[school] = average_squared_degree * (1 / average_degrees[school])
end

# friendship paradox ratio
paradox_ratios :: Dict{String, Float64} = Dict()
for (school, average_adjacent_degree) in average_adjacent_degrees
  paradox_ratios[school] = average_adjacent_degree / average_degrees[school]
end

# print(paradox_ratio)

## scatter plot where x = mean degree
## y = paradox ratio

if(keys(paradox_ratios) != 100)
  println("PROBLEMNMMM")
end

if(keys(average_degrees) != 100)
  println("PROBLEMNMMM")
end
print(length(keys(average_degrees)))
print(length(keys(paradox_ratios)))

println("+++++")
# average friendship paradox
summm = 0
for (k, para) in paradox_ratios
  global summm
  summm += para
end
print(summm / 100)

using Plots
gr()
plot(xlabel = "Mean degree ⟨kᵤ⟩", ylabel = "Friendship Paradox Ratio ⟨kᵥ⟩ / ⟨kᵤ⟩",
     title = "Friendship Paradox Ratio vs Mean Degree", legend = true, grid = true 
    )

for (school, paradox_ratio) in paradox_ratios
  average_degree = average_degrees[school]
  if (school == "Colgate88")
    scatter!([average_degree], [paradox_ratio], label = "Colgate", color = :red, markersize = 5)
  elseif (school == "Reed98")
    scatter!([average_degree], [paradox_ratio], label = "Reed", color = :orange, markersize = 5)
  elseif (school == "Mississippi66")
    scatter!([average_degree], [paradox_ratio], label = "Mississippi", color = :purple, markersize = 5)
  elseif (school == "Virginia63")
    scatter!([average_degree], [paradox_ratio], label = "Virginia", color = :teal, markersize = 5)
  elseif (school == "Berkeley13")
    scatter!([average_degree], [paradox_ratio], label = "UC Berkeley", color = :yellow, markersize = 5)

  else
    scatter!([average_degree], [paradox_ratio], label = "", color = :blue, markersize = 5)

  end
end 

hline!([1], label = "No Paradox", linestyle = :solid, color = :green)

savefig(current(), "ctest.pdf")







# <k_v> / <k_u> for each school where 
# <k_u> = 2m/n
# <k_v> = average_squared_degrees * (1 / <k_u>)
# average_squared_degrees = squared_degrees / n












#########




### c ###


#########

### d ###


#########


### e ###

# extra credit

########