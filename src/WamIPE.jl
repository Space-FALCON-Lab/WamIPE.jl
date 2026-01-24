module WamIPE

include("internal/Internal.jl")
include("Density.jl")
include("Live.jl")

export Density, Live, mean_density_profile, plot_global_mean_profile, 
plot_global_mean_profile_makie

end # module WamIPE
