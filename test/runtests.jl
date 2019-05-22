using Deldir
using Test

import DataFrames


@testset "A Voronoi cell for each generator" begin
    N = rand(5:15)
    x = rand(N)
    y = rand(N)
    summ = deldir(x, y)[3]

    @test DataFrames.nrow(summ) == N
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

    A = voronoiarea(x, y, [0., 2., -0.5, 1.])
    @test sum(A) ≈ 3 atol = 0.001
end


@testset "Error when points are outside window" begin
    x = [-0.1, 0.5]
    y = rand(2)

    @test_throws DomainError deldir(x, y)
end
