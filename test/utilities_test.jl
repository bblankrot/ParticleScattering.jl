@testset "utilities" begin
    s = rounded_star(1, 0.1, 7, 100)
    @test ParticleScattering.pInPolygon([0.0 0.0], s.ft)
    @test !ParticleScattering.pInPolygon([0.84 0.35], s.ft)
    @test ParticleScattering.pInPolygon([0.84 0.35], s.ft+0.1)

    sp = ScatteringProblem([s],[1],[0.1 0.1], [15π])
    border = find_border(sp)
    border2 = find_border(sp,[0.1 0.1])
    @test border == border2
    border3 = find_border(sp,[0.1 0.9999])
    @test border != border3
end
