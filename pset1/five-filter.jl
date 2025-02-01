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



names :: Vector{String} = get_names()

# collect done so far names
dones :: Set{String} = Set()


dir_path_data = "./data/"



for file in readdir(dir_path_data)
  push!(dones, file)
end 

names = filter(n -> !(n in dones), names)

println(length(names))



println(names)



# output as a list


