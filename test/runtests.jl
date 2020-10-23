using Deldir
using Test

import DataFrames


@testset "Deldir output are expected dataframes" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    del, vor, summ = deldir(x, y)

    @test DataFrames.names(del) == ["x1", "y1", "x2", "y2", "ind1", "ind2"]

    @test DataFrames.names(vor) == ["x1", "y1", "x2", "y2", "ind1", "ind2", "bp1", "bp2", "thirdv1", "thirdv2"]

    @test DataFrames.names(summ) == ["x", "y", "ntri", "del_area", "n_tside", "nbpt", "vor_area"]
    @test DataFrames.nrow(summ) == N

    @test summ[!, :x] == x
    @test summ[!, :y] == y
end


@testset "Delaunay triangle corners are indexed correctly" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    del = deldir(x, y)[1]

    @test del[!, :x1] == x[del[!, :ind1]]
    @test del[!, :x2] == x[del[!, :ind2]]

    @test del[!, :y1] == y[del[!, :ind1]]
    @test del[!, :y2] == y[del[!, :ind2]]
end


@testset "Line segments of Voronoi cells are within window" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    vor = deldir(x, y)[2]

    @test all(0 .<= vor[!, :x1] .<= 1)
    @test all(0 .<= vor[!, :y1] .<= 1)
    @test all(0 .<= vor[!, :x2] .<= 1)
    @test all(0 .<= vor[!, :y2] .<= 1)
end


@testset "Area of Voronoi cells sum to area of window" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    A = voronoiarea(x, y)

    @test sum(A) ≈ 1 atol = 0.001

    rw = [-rand(); 1 + rand(); -rand(); 1 + rand()]
    rw_area = (rw[2] - rw[1])*(rw[4] - rw[3])

    A = voronoiarea(x, y, rw)

    @test sum(A) ≈ rw_area atol = 0.001
end


@testset "Errors with inappropriate input" begin
    @testset "Error when points are outside window" begin
        x = [-rand(), rand()]
        y = rand(2)
    
        @test_throws DomainError deldir(x, y)
    end
    
    @testset "Error when number of x's and y's are not equal" begin
        x = rand(rand(2:7))
        y = rand(rand(8:12))
    
        @test_throws DimensionMismatch deldir(x, y)
    end
end


@testset "Handle Fortran errors" begin
    # Data from GitHub issue #17
    x = [0.4, 0.3, 0.5, 0.2406, 0.2964, 0.5498, 0.2332, 0.3, 0.5041, 0.0824, 0.0594, 0.0126, 0.4385, 0.3575, 0.7737, 0.1, 0.1997, 0.6806, 0.8219, 0.0098, 0.4568, 0.0136]
    y = [0.3856, 0.5588, 0.0, 0.0725, 0.0433, 0.0025, 0.0771, 0.2124, 0.0, 0.2251, 0.7363, 0.3885, 0.0038, 0.0207, 0.0816, 0.2124, 0.1002, 0.0338, 0.3856, 0.4017, 0.0019, 0.616]

    @test_logs (:info, "Fortran error 4. Increasing madj to 24") voronoiarea(x, y)
end

