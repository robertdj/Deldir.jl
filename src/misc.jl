"""
	remove_duplicates(x::Vector, y::Vector)

Remove duplicate tuples `(x[i],y[i])` from the vectors `x` and `y`.
"""
function remove_duplicates(x::Vector, y::Vector)
	points = hcat(x, y)
	unique_points = unique(points, dims = 1)

	return unique_points[:, 1], unique_points[:, 2]
end


"""
	voronoiarea(x::Vector, y::Vector, rw) -> Vector

Compute the area of each Voronoi cell of the generators `(x[i],y[i])` in the vectors `x` and `y`.

`rw` is the boundary window.
"""
function voronoiarea(x::Vector, y::Vector, rw::Vector=[0.0; 1.0; 0.0; 1.0])
    da = DeldirArguments(x, y, rw, 1e-9)
	deldirwrapper!(da)

    npd = Int64(da.npd[1])
	da.dirsum[2*npd + 1:3*npd]
end

"""
	edges(D) -> Vector, Vector

Collect the edges of a dataframe in vectors that are ready to be plotted.
"""
function edges(D::DataFrame)
	x1 = D[:x1]
	y1 = D[:y1]
	x2 = D[:x2]
	y2 = D[:y2]

	N = size(D, 1)
	x = Array{Float64}(undef, 3*N)
	y = similar(x)
	
	nx = 0
	for n = 1:N
		x[nx+=1] = x1[n]
		y[nx] = y1[n]

		x[nx+=1] = x2[n]
		y[nx] = y2[n]

		x[nx+=1] = NaN
		y[nx] = NaN
	end

	return x, y
end

