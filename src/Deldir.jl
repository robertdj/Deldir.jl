module Deldir

using DataFrames
using Deldir_jll

export
	deldir,
	voronoiarea,
	edges

include("wrapper.jl")
include("misc.jl")

end # module
