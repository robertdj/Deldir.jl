using Deldir
using Test

import DataFrames


@testset "Deldir output are expected dataframes" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    del, vor, summ = deldir(x, y)

    @test DataFrames.names(del) == [:x1, :y1, :x2, :y2, :ind1, :ind2]

    @test DataFrames.names(vor) == [:x1, :y1, :x2, :y2, :ind1, :ind2, :bp1, :bp2]

    @test DataFrames.names(summ) == [:x, :y, :ntri, :del_area, :n_tside, :nbpt, :vor_area]
    @test DataFrames.nrow(summ) == N

    @test summ[:x] == x
    @test summ[:y] == y
end


@testset "Delaunay triangle corners are indexed correctly" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    del = deldir(x, y)[1]

    @test del[:x1] == x[del[:ind1]]
    @test del[:x2] == x[del[:ind2]]

    @test del[:y1] == y[del[:ind1]]
    @test del[:y2] == y[del[:ind2]]
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
    
    @testset "Error number of x's and y's are not equal" begin
        x = rand(rand(2:7))
        y = rand(rand(8:12))
    
        @test_throws DimensionMismatch deldir(x, y)
    end
end
