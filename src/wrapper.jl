mutable struct DeldirArguments
    # The variables for the master Fortran subroutine
    # The variable names are copied from the R package
    x::Vector{Float64}
    y::Vector{Float64}
    rw::Vector{Float64}
    npd::Vector{Int32}
    ntot::Vector{Int32}
    nadj::Vector{Int32}
    madj::Vector{Int32}
    tx::Vector{Float64}
    ty::Vector{Float64}
    epsilon::Float64
    delsgs::Vector{Float64}
    ndel::Vector{Int32}
    delsum::Vector{Float64}
    dirsgs::Vector{Float64}
    ndir::Vector{Int32}
    dirsum::Vector{Float64}
    nerror::Vector{Int32}

    # The variables for sorting the points
    indices::Vector{Int32}
    reverse_indices::Vector{Int32}

    function DeldirArguments(x, y, rw, npd, ntot, nadj, madj, tx, ty, epsilon, delsgs, ndel, 
                             delsum, dirsgs, ndir, dirsum, nerror, indices, reverse_indices)
        if length(x) != length(y)
            throw(DimensionMismatch("Coordinate vectors must be of equal length"))
        end

        if epsilon < eps(Float64)
            throw(DomainError(epsilon, "Must be at least `eps(Float64)`"))
        end

        min_x, max_x = extrema(x)
        min_y, max_y = extrema(y)
        if min_x < rw[1] || max_x > rw[2] || min_y < rw[3] || max_y > rw[4]
            throw(DomainError(rw, "Boundary window is too small"))
        end

        new(x, y, rw, npd, ntot, nadj, madj, tx, ty, epsilon, delsgs,
            ndel, delsum, dirsgs, ndir, dirsum, nerror, indices, reverse_indices)
    end
end


function DeldirArguments(x, y, rw, epsilon)
    unique_x, unique_y = remove_duplicates(x, y)
    num_points = length(unique_x)

    # Dummy points: Ignored!
    ndm = 0
    npd = num_points + ndm

    # The total number of points
    ntot = npd + 4

    indices, reverse_indices = sortperm_points!(unique_x, unique_y, rw)
    X = [zeros(4); unique_x; zeros(4)]
    Y = [zeros(4); unique_y; zeros(4)]

    # Set up fixed dimensioning constants
    ntdel = 4*npd
    ntdir = 3*npd

    # Set up dimensioning constants which might need to be increased
    madj_val = max(20, ceil(Int32, 3*sqrt(ntot)))
    madj = Int32[madj_val]
    tadj = (madj_val + 1)*(ntot + 4)
    ndel = Int32[madj_val*(madj_val + 1)/2]
    tdel = 6*ndel[]
    ndir = copy(ndel)
    tdir = 10*ndir[]

    nadj   = zeros(Int32, tadj)
    tx     = zeros(Float64, npd)
    ty     = zeros(Float64, npd)
    delsgs = zeros(Float64, tdel)
    delsum = zeros(Float64, ntdel)
    dirsgs = zeros(Float64, tdir)
    dirsum = zeros(Float64, ntdir)
    nerror = Int32[1]

    DeldirArguments(X, Y, rw, [Int32(npd)], [Int32(ntot)], nadj, madj, tx, ty, epsilon, 
                    delsgs, ndel, delsum, dirsgs, ndir, dirsum, nerror, indices, reverse_indices)
end


function sortperm_points!(x, y, rw)
    n = length(x)
    if n != length(y)
        DimensionMismatch("x and y must have the same length")
    end

    tx = similar(x)
    ty = similar(y)

    indices = Vector{Int32}(undef, n)
    reverse_indices = similar(indices)
    ilst = similar(indices)
    nerror = Int32[1]

    ccall((:binsrt_, Deldir_jll.libdeldir), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, 
           Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}),
          x, y, rw, n, indices,
          reverse_indices, tx, ty, ilst, nerror
    )

    if nerror[] > 0
        error("Mismatch between number of points and number of sorted points")
    end

    return indices, reverse_indices
end


function error_handling!(da::DeldirArguments)
    error_number = da.nerror[]

    if error_number == 4
        madj_val = ceil(Int32, 1.2*da.madj[])
        da.madj = Int32[madj_val]

        tadj = (madj_val + 1)*(da.ntot[] + 4)
        resize!(da.nadj, tadj)

        ndel_val = max(da.ndel[], div(madj_val*(madj_val + 1), 2))
        da.ndel = Int32[ndel_val]
        da.ndir = copy(da.ndel)

        tdel = 6*ndel_val
        resize!(da.delsgs, tdel)

        tdir = 8*ndel_val
        resize!(da.dirsgs, tdir)

        @info "Fortran error $(error_number). Increasing madj to $(da.madj[])"
    elseif error_number == 12
        error("Vertices of triangle are collinear")
    elseif error_number == 14 || error_number == 15
        ndel_val = ceil(Int32, 1.2*ndel[])
        da.ndel = Int32[ndel_val]
        da.ndir = copy(ndel)

        tdel = 6*ndel_val
        resize!(da.delsgs, tdel)

        tdir = 8*ndel_val
        resize!(da.dirsgs, tdir)

        @info "Fortran error $(error_number). Increasing ndel & ndir to $(da.ndel[])"
    elseif error_number > 1
        error("From `deldir` Fortran, nerror = ", error_number)
    end

    return(da)
end


function get_delaunay(da::DeldirArguments)
    num_del = Int64(da.ndel[])
    delsgs  = reshape(da.delsgs[1:6*num_del], 6, num_del) |> transpose

    ind1 = round.(Int, delsgs[:, 5])
    ind2 = round.(Int, delsgs[:, 6])

    DataFrames.DataFrame(
        x1 = delsgs[:, 1],
        y1 = delsgs[:, 2],
        x2 = delsgs[:, 3],
        y2 = delsgs[:, 4],
        ind1 = da.reverse_indices[ind1],
        ind2 = da.reverse_indices[ind2]
    )
end


function get_voronoi(da::DeldirArguments)
    num_dir = Int64(da.ndir[])
    vor  = reshape(da.dirsgs[1:10*num_dir], 10, num_dir) |> transpose

    ind1 = round.(Int, vor[:, 5])
    ind1 = da.reverse_indices[ind1]

    ind2 = round.(Int, vor[:, 6])
    ind2 = da.reverse_indices[ind2]

    thirdv1 = Int.(vor[:, 9])
    idx1 = thirdv1 .>= 0
    thirdv1[idx1] = da.reverse_indices[thirdv1[idx1]]

    thirdv2 = Int.(vor[:, 10])
    idx2 = thirdv2 .>= 0
    thirdv2[idx2] = da.reverse_indices[thirdv2[idx2]]

    DataFrames.DataFrame(
        x1 = vor[:, 1],
        y1 = vor[:, 2],
        x2 = vor[:, 3],
        y2 = vor[:, 4],
        ind1 = ind1,
        ind2 = ind2,
        bp1 = vor[:, 7] .== 1,
        bp2 = vor[:, 8] .== 1,
        thirdv1 = thirdv1,
        thirdv2 = thirdv2
    )
end


function get_summary(da::DeldirArguments)
    npd = Int64(da.npd[1])

    delsum = reshape(da.delsum, npd, 4)
    delsum = delsum[da.indices, :]

    dirsum = reshape(da.dirsum, npd, 3)
    dirsum = dirsum[da.indices, :]
    
    DataFrames.DataFrame(
        x = delsum[:, 1],
        y = delsum[:, 2],
        ntri = round.(Int, delsum[:, 3]),
        del_area = delsum[:, 4],
        n_tside = round.(Int, dirsum[:, 1]),
        nbpt = round.(Int, dirsum[:, 2]),
        vor_area = dirsum[:, 3],
    )
end


"""
deldirwrapper!(da::DeldirArguments)

Wrapper for the Fortran code that returns the output undigested.
"""
function deldirwrapper!(da::DeldirArguments)
    # Call Fortran routine
    while da.nerror[] >= 1
        ccall((:master_, Deldir_jll.libdeldir), Cvoid,
            (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, 
             Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64},
             Ref{Float64}, Ref{Float64}, Ref{Int32}, 
             Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Float64}, Ref{Int32}),
            da.x, da.y, da.rw, da.npd, da.ntot,
            da.nadj, da.madj, da.tx, da.ty,
            da.epsilon, da.delsgs, da.ndel,
            da.delsum, da.dirsgs, da.ndir, da.dirsum, da.nerror
        )

        error_handling!(da)
    end

    return da
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
- The `thirdv1` and `thirdv2` columns are the indices of the respective third vertices of the Delaunay triangle whose circumcentres constitute the corresponding endpoints of the edge under consideration.

###### `summary`

- The `x` and `y` entries of each row are the coordinates of the points in the set being triangulated.
- The `ntri` entry are the number of Delaunay triangles emanating from the point.
- The `del_area` entry is `1/3` of the total area of all the Delaunay triangles emanating from the point.
- The `n_tside` entry is the number of sides — within the rectangular window — of the Voronoi cell surrounding the point.
- The `nbpt` entry is the number of points in which the Voronoi cell intersects the boundary of the rectangular window.
- The `vor_area` entry is the area of the Voronoi cell surrounding the point.
"""
function deldir(x::Vector{Float64}, y::Vector{Float64}, rw::Vector = [0.0; 1.0; 0.0; 1.0], epsilon = 1e-9)
    da = DeldirArguments(x, y, rw, epsilon)
    deldirwrapper!(da)

    del = get_delaunay(da)
    vor = get_voronoi(da)
    summary = get_summary(da)

    return del, vor, summary
end
