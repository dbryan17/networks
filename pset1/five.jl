

# import 
dir_path = "./facebook100txt/"

get_names = () -> 
  # get names of ones that acutally contain data
  filter(file -> endswith(file, ".txt") && 
                          !endswith(file, "_attr.txt") && 
                          !occursin("readme", file),
                 readdir(dir_path))



names = get_names()

# is 100 which is correct
print(length(names))

edge_set = Set()

# check if edges are counted twice
open(dir_path * names[1]) do f 

  for line in eachline(f)
    nodes = split(line)
    if(length(nodes) != 2)
      print("errr")
    end
    edge = tuple(min(nodes[1], nodes[2]), max(nodes[1], nodes[2]))
    if edge in edge_set
      print("DUPLICATED")
    else
      push!(edge_set, edge)
    end

  end
end








### a ###





#########


### b ###


### c ###


#########

### d ###


#########


### e ###

# extra credit

########