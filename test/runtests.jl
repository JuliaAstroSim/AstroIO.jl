using Test, Unitful, UnitfulAstro
using FileIO

using PhysicalParticles

using AstroIO

header, data = read_gadget2("gassphere_littleendian.gadget2") # 1472 gas particles

@testset "Gadget" begin
    result = print(header)
    @test isnothing(result)

    @test length(data) == 1472

    @test write_gadget2("testGadget.gadget2", header, data)

    @test write_gadget2("testGadgetHeaderGeneration.gadget2", data)

    # format2
    @test write_gadget2_format2("gadget2.format2", header, data)

    h, d = read_gadget2("gadget2.format2")
    @test h.npart[1] == 1472

    pos = read_gadget2_pos("gadget2.format2")
    @test length(pos) == 1472

    pos = read_gadget2_pos("gassphere_littleendian.gadget2")
    @test length(pos) == 1472

    h, d = read_gadget2("pot_acc.format2.gadget2", acc = true, pot = true)
    write_gadget2_format2("pot_acc.format2.test.gadget2", h, d, acc = true, pot = true)
    h, d = read_gadget2("pot_acc.format2.test.gadget2", acc = true, pot = true)
    @test !iszero(norm(average(d, :Acc)))
    @test !iszero(norm(average(d, :Potential)))
end

@testset "FileIO" begin
    h, d = load("gassphere_littleendian.gadget2")
    @test length(d) == 1472

    @test isnothing(save("testFileIO.gadget2", h, d))
end

@testset "CSV" begin
    stars2d = [Star2D() for i = 1:10]
    @test write_csv("testcsvStar2D", stars2d, nothing)

    stars = [Star(uAstro) for i = 1:10]
    @test write_csv("testcsvStar", stars, uAstro)

    data = [[Star() for i = 1:10]; [Ball() for i = 1:10]]
    @test write_csv("testcsvGeneral", data, nothing)

    @test write_ramses("ramses.csv", stars)
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

@testset "Tools" begin
    f = open("test_renamesuffixs.test", "w")
    close(f)
    
    renamesuffixs(pwd(), "test_rename", ".ok")

    @test isfile("test_renamesuffixs.ok")
end

#=
include("AstroIO.jl\\src\\AstroIO.jl"); using .AstroIO

h, d = read_gadget2("AstroIO.jl\\test\\gassphere_littleendian.gadget2")
write_gadget2("AstroIO.jl\\test\\testGadget.gadget2", h, d)

h, d = load("AstroIO.jl\\test\\gassphere_littleendian.gadget2")
save("AstroIO.jl\\test\\testFileIO.gadget2", h, d)
=#