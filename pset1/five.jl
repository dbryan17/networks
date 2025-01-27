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



function count_dups(filename) 
  # check if edges are counted twice
  open(dir_path * filename) do f 
    dup_count = 0
    total_count = 0
    edge_set = Set()
    global total_count
    global dup_count
    global edge_set
    for line in eachline(f)
      nodes = split(line)
      if(length(nodes) != 2)
        print("errr")
      end
      total_count += 1
      node1, node2 = parse(Int, nodes[1]), parse(Int, nodes[2])
      edge = tuple(min(node1, node2), max(node1, node2))
      if edge in edge_set
        dup_count += 1
      else
        push!(edge_set, edge)
      end
    end
    if (total_count / 2 != dup_count)
      print("duplicates are not half!!")
    end

  end
 
end  


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





# ok, so all undirected, and data is wel formed
function check_undirected()
  for name in get_names()
    count_dups(name)
  end
end


# check_undirected()


### a ###


edges :: Dict{String, Set{Tuple{Int, Int}}} = get_all_edges()

edges_counts :: Dict{String, Int} = Dict()
for (k,v) in edges
  edges_counts[k] = length(v)
end

total_nodes :: Dict{String, Int} = Dict()
for (k,v) in edges
  nodes = Set()
  for (n1, n2) in v
    push!(nodes, n1)
    push!(nodes, n2)
  end
  total_nodes[k] = length(nodes)
end

avgs = []
for (k,v) in total_nodes
  # (2 * |e|) / |nodes|
  avg = (2 * edges_counts[k]) / v
  push!(avgs, avg)

end


if length(avgs) != 100
  print("probem")

end


using Plots
gr()
p = histogram(avgs, binwidth=5, title="Average Degree Distribution in Facebook 100", xlabel="Average Degree", ylabel="Frequency", legend=false)
savefig(p, "aa.pdf")


#########


### b ###




### c ###


#########

### d ###


#########


### e ###

# extra credit

########