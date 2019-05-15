using Deldir
using Test

@testset "Area of Voronoi cells sum to area of window" begin
    N = 10
    x = rand(N)
    y = rand(N)

    A = voronoiarea(x, y)
    @test sum(A) ≈ 1
end

