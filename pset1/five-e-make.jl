figure_max_total_size :: Dict{String, Tuple{Int, Int}} = Dict()
figure_avg_comp_size :: Dict{String, Tuple{Float64, Int}} = Dict()

dir_path = "./data/"

for name in readdir(dir_path)
  open(dir_path * name) do f 
    short_name :: String = ""
    avg_dist :: Float64 = 0
    total_nodes :: Int = 0
    largest_comp :: Int = 0
    max_dist :: Int = 0
    for (line_number, line) in enumerate(eachline(f))
      if line_number == 1
        short_name = line

      end
      if line_number == 3 
        avg_dist = parse(Float64, line)
      end

      if line_number == 5
        largest_comp = parse(Int, line)
      end

      if line_number == 7 
        max_dist = parse(Int, line)
      end 


      if line_number == 9 
        total_nodes = parse(Int, line)
      end

    end 

    figure_avg_comp_size[short_name] = (avg_dist, largest_comp)

    figure_max_total_size[short_name] = (max_dist, total_nodes)
  
  
    
  end


end

max_fig_x = []
max_fig_y = []
diameters = []
for (name, (max, total)) in figure_max_total_size
  push!(max_fig_x, total)
  push!(max_fig_y, max)
  push!(diameters, max)
end


avg_fig_x = []
avg_fig_y = []

for (name, (avg, largest)) in figure_avg_comp_size
  push!(avg_fig_x, largest)
  push!(avg_fig_y, avg)
end


using Plots
gr()
max_fig = scatter(max_fig_x, max_fig_y, title="Diameter as a function of total network size", xlabel="network size", ylabel="diameter", color=:blue, markersize=5, legend=false)
max_fig_hist = histogram(diameters, binwidth=1, title="Diameter Distrubtion", xlabel="Diameter", ylabel="Frequency", legend=false)
avg_fig = scatter(avg_fig_x, avg_fig_y, title="Mean geodesic distance vs largest component size", xlabel="largest component size", ylabel="mean geodesic distance", color=:blue, markersize=5, legend=false)
savefig(max_fig, "e-maxfig-new.pdf")
savefig(avg_fig, "e-avgfig-new.pdf")
savefig(max_fig_hist, "e-maxfig-hist.pdf")
