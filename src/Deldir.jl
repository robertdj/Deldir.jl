module Deldir

using DataFrames

# Get path for libdeldir.so as libdeldir
depsfile = joinpath(@__DIR__, "..", "deps", "deps.jl")
if isfile(depsfile)
	include(depsfile)
end

export
	deldir,
	voronoiarea,
	edges

include("wrapper.jl")
include("misc.jl")

end # module
