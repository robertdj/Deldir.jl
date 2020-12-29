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

Compute the area of each Voronoi cell of the generators `(x[i], y[i])` in the vectors `x` and `y`.

`rw` is the boundary window.
"""
function voronoiarea(x::Vector, y::Vector, rw::Vector=[0.0; 1.0; 0.0; 1.0])
    da = DeldirArguments(x, y, rw, 1e-9)
	deldirwrapper!(da)

    npd = Int64(da.npd[1])
	unordered_area = @view da.dirsum[2*npd + 1:3*npd]

	area = unordered_area[da.indices]
	return area
end


"""
	edges(D) -> Vector, Vector

Collect the edges of a dataframe in vectors that are ready to be plotted.
"""
function edges(D::DataFrame)
	x1 = D[!, "x1"]
	y1 = D[!, "y1"]
	x2 = D[!, "x2"]
	y2 = D[!, "y2"]

	N = DataFrames.nrow(D)
	x = Array{Float64}(undef, 3*N)
	y = similar(x)
	
	nx = 0
	for n in 1:N
        nx += 1
		x[nx] = x1[n]
		y[nx] = y1[n]

        nx += 1
		x[nx] = x2[n]
		y[nx] = y2[n]

        nx += 1
		x[nx] = NaN
		y[nx] = NaN
	end

	return x, y
end

