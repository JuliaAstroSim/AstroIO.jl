using Test, Unitful, UnitfulAstro
using FileIO

using PhysicalParticles

using AstroIO

header, data = read_gadget2("gassphere_littleendian.g2")

@testset "Gadget" begin
    @test length(data) == 1472

    @test write_gadget2("testGadget.g2", header, data)
end

@testset "FileIO" begin
    h, d = load("gassphere_littleendian.g2")
    @test length(d) == 1472

    @test save("testFileIO.g2", h, d)
end

@testset "CSV" begin
    stars2d = [Star2D() for i = 1:10]
    @test write_csv("testcsvStar2D", stars2d)

    stars = [Star() for i = 1:10]
    @test write_csv("testcsvStar", stars)

    gases2d = [SPHGas2D() for i=1:10]
    @test write_csv("testcsvSPHGas2D", gases2d)

    gases = [SPHGas() for i=1:10]
    @test write_csv("testcsvSPHGas", gases)

    data = [[Star() for i = 1:10]; [SPHGas() for i = 1:10]]
    @test write_csv("testcsvGeneral", data)
end

@testset "JLD2" begin
    d = [Star2D() for i = 1:10]
    @test write_gadget2_jld("testjldGadget.jld2", header, d)

    h, d = read_gadget2_jld("testjldGadget.jld2")
    @test length(d) == 10

    @test write_jld("testjldGeneral.jld2", d)

    d = read_jld("testjldGadget.jld2")
    @test length(d) == 10
end


#=
include("AstroIO.jl\\src\\AstroIO.jl"); using .AstroIO

h, d = read_gadget2("AstroIO.jl\\test\\gassphere_littleendian.g2")
write_gadget2("AstroIO.jl\\test\\testGadget.g2", h, d)

h, d = load("AstroIO.jl\\test\\gassphere_littleendian.g2")
save("AstroIO.jl\\test\\testFileIO.g2", h, d)
=#