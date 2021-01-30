import VoronoiCells
import Deldir
using GeometryBasics
using Plots

# Number of generators: Deldir does not go all the way
Nsmall = [1000; 2000:2000:30_000]
Nbig = 40_000:10_000:100_000
#= Nsmall = [1000; 2000:2000:10_000] =#
#= Nbig = [10_000, 20_000] =#
N = vcat(Nsmall, Nbig)

VCtime = Vector{Float64}(undef, length(N))
Dtime = similar(VCtime)
fill!(Dtime, NaN)

# Compile methods
Deldir.voronoiarea(rand(8), rand(8))
VoronoiCells.voronoicells(rand(8), rand(8), VoronoiCells.Rectangle(Point2(0, 0), Point2(1, 1)))

# Simulation window
W = [0.0 ; 10.0 ; 0.0 ; 10.0]
rect = VoronoiCells.Rectangle(Point2(0, 0), Point2(10, 10))

for (idx, n) in enumerate(N)
	println(n)

	local x = W[1] .+ W[2]*rand(n)
	local y = W[3] .+ W[4]*rand(n)

    if n <= Nsmall[end]
		Dtime[idx] = @elapsed Deldir.voronoiarea(x, y, W)
	end

    local points = [Point2(x[k], y[k]) for k in 1:n]

    VCtime[idx] = @elapsed begin
        tess = VoronoiCells.voronoicells(points, rect)
        VoronoiCells.voronoiarea(tess)
    end
end


# ------------------------------------------------------------

x = div.(N, 1000)
scatter(x, Dtime, label = "Deldir")
scatter!(x, VCtime, label = "VoronoiCells")
plot!(xlabel = "number of points in 1000s", ylabel = "time in seconds", xticks = 0:20:100, yticks = 0:5:20)
plot!(tickfont = font(13), legendfont = font(12))

savefig("comparison.png")

