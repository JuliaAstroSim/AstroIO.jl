using Test, Unitful, UnitfulAstro

using PhysicalParticles

using AstroIO

@testset "Gadget" begin
    h, d = read_gadget2("gassphere_littleendian.dat")
    @test length(d.gases) == 1472
end