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


check_undirected()