using Test, Unitful, UnitfulAstro

using PhysicalParticles

using AstroIO

@testset "Gadget" begin
    h, d = read_gadget2("gassphere_littleendian.g2")
    @test length(d.gases) == 1472

    @test write_gadget2("test.g2", h, d)
end

@testset "CSV" begin
    d = Dict("Stars" => [Star() for i = 1:10], "SPHGases" => [SPHGas() for i=1:10])
    @test write_csv("csvtest", d)
    @test write_csv("csvtest", d, seperate = true)
end