using Test, Unitful, UnitfulAstro

using PhysicalParticles

using AstroIO

@testset "Gadget" begin
    h, d = read_gadget2("gassphere_littleendian.g2")
    @test length(d.gases) == 1472

    @test write_gadget2("test.g2", h, d)
end