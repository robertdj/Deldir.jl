# Disclaimer: The macros initialize, error_handling and finalize use variable 
# names that are copied verbatim from the R package. 


"""
	remove_duplicates(x::Vector, y::Vector)

Remove duplicate tuples `(x[i],y[i])` from the vectors `x` and `y`.
"""
function remove_duplicates(x::Vector, y::Vector)
	points = [x y]
	unique_points = unique(points, dims = 1)

	return unique_points[:, 1], unique_points[:, 2]
end

"""
	initialize()

Set up input for the deldir Fortran routine
"""
macro initialize()
	esc(quote
		# According to the documentation in the R package: 
		# "'sort' would get used only in a de-bugging process"
		# Therefore it is not an argument to the wrapper
		sort = 1

		x, y = remove_duplicates(x, y)
		num_points = length(x)

		# Dummy points: Ignored!
		ndm = 0
		npd = num_points + ndm

		# The total number of points
		ntot = npd + 4

		X = [zeros(4); x; zeros(4)]
		Y = [zeros(4); y; zeros(4)]

		# Set up fixed dimensioning constants
		ntdel = 4*npd
		ntdir = 3*npd

		# Set up dimensioning constants which might need to be increased
		madj = max(20, ceil(Int32, 3*sqrt(ntot)))
		tadj = (madj + 1)*(ntot + 4)
		ndel = Int32[madj*(madj + 1)/2]
		tdel = 6*ndel[]
		ndir = copy(ndel)
		tdir = 8*ndir[]

		nadj   = zeros(Int32, tadj)
        tx     = zeros(Float64, npd)
        ty     = zeros(Float64, npd)
		delsgs = zeros(Float64, tdel)
		delsum = zeros(Float64, ntdel)
		dirsgs = zeros(Float64, tdir)
		dirsum = zeros(Float64, ntdir)
		nerror = Int32[1]

		@allocate
	end)
end

"""
	allocate()

Allocate input to be modified by the deldir Fortran routine
"""
macro allocate()
	esc(quote
		nadj   = zeros(Int32, tadj)
        tx     = zeros(Float64, npd)
        ty     = zeros(Float64, npd)
		delsgs = zeros(Float64, tdel)
		delsum = zeros(Float64, ntdel)
		dirsgs = zeros(Float64, tdir)
		dirsum = zeros(Float64, ntdir)
		nerror = Int32[1]
	end)
end

"""
	error_handling()

Handle errors from the deldir Fortran routine
"""
macro error_handling()
	esc(quote
		if nerror[] == 4
			madj = ceil(Int32, 1.2*madj)
			tadj = (madj + 1)*(ntot + 4)
			ndel = max(ndel, div(madj*(madj + 1), 2))
			tdel = 6*ndel[]
			ndir = copy(ndel)
			tdir = 8*ndir[]

			@allocate
		elseif nerror[] == 14 || nerror[] == 15
			ndel = ceil(Int32, 1.2*ndel)
			tdel = 6*ndel[]
			ndir = copy(ndel)
			tdir = 8*ndir[]

			@allocate
		elseif nerror[] > 1
			error("From `deldir` Fortran, nerror = ", nerror[])
		end
	end)
end

"""
	finalize()

Process output from the deldir Fortran routine
"""
macro finalize()
	esc(quote
		num_del = Int64(ndel[])
        delsgs  = transpose(reshape(delsgs[1:6*num_del], 6, num_del))
		num_dir = Int64(ndir[])
        dirsgs  = transpose(reshape(dirsgs[1:8*num_dir], 8, num_dir))
		delsum  = reshape(delsum, npd, 4)
		dirsum  = reshape(dirsum, npd, 3)
		allsum  = hcat(delsum, dirsum)
	end)
end

"""
	deldirwrapper(x::Vector{Float64}, y::Vector{Float64}; ...)

Wrapper for the Fortran code that returns the output rather undigested.
"""
function deldirwrapper(x::Vector{Float64}, y::Vector{Float64}, 
                       rw::Vector = [0.0; 1.0; 0.0; 1.0]; epsilon::Float64 = 1e-9)

	if length(x) != length(y)
        throw(DimensionMismatch("Coordinate vectors must be of equal length"))
    end

    if epsilon < eps(Float64)
        throw(DomainError())
    end

	if minimum(x) < rw[1] || maximum(x) > rw[2] && minimum(y) < rw[3] && maximum(y) > rw[4] 
        throw(DomainError("Boundary window is too small"))
    end

	@initialize

	# Call Fortran routine
	while nerror[] >= 1
		ccall((:master_, Deldir_jll.libdeldir), Cvoid,
		    (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, 
             Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
             Ref{Float64}, Ref{Float64}, Ref{Int32}, 
             Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Int32}),
		    X, Y, float(rw), npd, ntot, 
            nadj, madj, tx, ty, 
            epsilon, delsgs, ndel, 
            delsum, dirsgs, ndir, dirsum, nerror
		)

		@error_handling
	end

	@finalize

	return delsgs, dirsgs, allsum
end

"""
	deldir(x::Vector, y::Vector; ...)

Compute the Delaunay triangulation and Voronoi tesselation of the 2D points with x-coordinates `x` and y-coordinates `y`.

Optional arguments are

- `rw`: Boundary rectangle specified as a vector with `[xmin, xmax, ymin, ymax]`. By default, `rw` is the unit rectangle.
- `epsilon`: A value of epsilon used in testing whether a quantity is zeros, mainly in the context of whether points are collinear.
If anomalous errors arise, it is possible that these may averted by adjusting the value of `epsilon` upward or downward.
By default, `epsilon = 1e-9`.

The output are three `DataFrame`s:

###### `delsgs`

- The `x1`, `y1`, `x2` & `y2` entires are the coordinates of the points joined by an edge of a Delaunay triangle.
- The `ind1` and `ind2` entries are the indices of the two points which are joined.

###### `vorsgs`

- The `x1`, `y1`, `x2` & `y2` entires are the coordinates of the endpoints of one the edges of a Voronoi cell.
- The `ind1` and `ind2` entries are the indices of the two points, in the set being triangulated, which are separated by that edge
- The `bp1` entry indicates whether the first endpoint of the corresponding edge of a Voronoi cell is a boundary point (a point on the boundary of the rectangular window). 
Likewise for the `bp2` entry and the second endpoint of the edge.

###### `summary`

- The `x` and `y` entries of each row are the coordinates of the points in the set being triangulated.
- The `ntri` entry are the number of Delaunay triangles emanating from the point.
- The `del_area` entry is `1/3` of the total area of all the Delaunay triangles emanating from the point.
- The `n_tside` entry is the number of sides — within the rectangular window — of the Voronoi cell surrounding the point.
- The `nbpt` entry is the number of points in which the Voronoi cell intersects the boundary of the rectangular window.
- The `vor_area` entry is the area of the Voronoi cell surrounding the point.
"""
function deldir(x::Vector{Float64}, y::Vector{Float64}; args...)
	del, vor, summ = deldirwrapper(x, y; args...)

    del_df = DataFrames.DataFrame(
        [Float64, Float64, Float64, Float64, Int, Int], 
        [:x1, :y1, :x2, :y2, :ind1, :ind2], 
        size(del, 1)
    )
	del_df[!, :x1]   = del[:, 1]
	del_df[!, :y1]   = del[:, 2]
	del_df[!, :x2]   = del[:, 3]
	del_df[!, :y2]   = del[:, 4]
	del_df[!, :ind1] = round.(Int, del[:, 5])
	del_df[!, :ind2] = round.(Int, del[:, 6])

    vor_df = DataFrames.DataFrame(
        [Float64, Float64, Float64, Float64, Int, Int, Bool, Bool], 
        [:x1, :y1, :x2, :y2, :ind1, :ind2, :bp1, :bp2], 
        size(vor, 1)
    )
	vor_df[!, :x1]   = vor[:, 1]
	vor_df[!, :y1]   = vor[:, 2]
	vor_df[!, :x2]   = vor[:, 3]
	vor_df[!, :y2]   = vor[:, 4]
	vor_df[!, :ind1] = round.(Int, vor[:, 5])
	vor_df[!, :ind2] = round.(Int, vor[:, 6])
	vor_df[!, :bp1]  = vor[:, 7] .== 1
	vor_df[!, :bp2]  = vor[:, 8] .== 1

    summary_df = DataFrames.DataFrame(
        [Float64, Float64, Int, Float64, Int, Int, Float64], 
        [:x, :y, :ntri, :del_area, :n_tside, :nbpt, :vor_area],
        size(summ, 1)
    )
	summary_df[!, :x]        = summ[:, 1]
	summary_df[!, :y]        = summ[:, 2]
	summary_df[!, :ntri]     = round.(Int, summ[:, 3])
	summary_df[!, :del_area] = summ[:, 4]
	summary_df[!, :n_tside]  = round.(Int, summ[:, 5])
	summary_df[!, :nbpt]     = round.(Int, summ[:, 6])
	summary_df[!, :vor_area] = summ[:, 7]

    del_df, vor_df, summary_df
end

