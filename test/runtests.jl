using Test, Unitful, UnitfulAstro
using FileIO

using PhysicalParticles

using AstroIO

@testset "GadgetFormat1" begin
    h, d = read_gadget2("gassphere_littleendian.gadget2", uAstro) # 1472 gas particles

    @test length(d) == 1472
    @test sum(h.npart) == 1472
    @test d.Mass[1] == 6.7934785f6u"Msun"
    @test d.Pos[1] == PVector(-0.07133729f0, -0.35668644f0, -0.92738473f0, u"kpc")

    @test write_gadget2("testGadget.gadget2", h, d)

    @test write_gadget2("testGadgetHeaderGeneration.gadget2", d)

    pos = read_gadget2_pos("gassphere_littleendian.gadget2", uAstro, uGadget2)
    @test length(pos) == 1472
    @test pos[1] == PVector(-0.07133729010820389,
                            -0.35668644309043884,
                            -0.9273847341537476, u"kpc")


end

@testset "GadgetFormat2" begin
    h, d = read_gadget2("pot_acc.format2.gadget2", uAstro, uGadget2)
    @test h.npart[5] == 1000
    @test d.Mass[1] == 100u"Msun"

    h, d = read_gadget2("pot_acc.format2.gadget2", nothing, uGadget2)
    @test d.Mass[1] == 1.0f-8

    h, d = read_gadget2("pot_acc.format2.gadget2", nothing, nothing)
    @test d.Mass[1] == 1.0f-8

    # getindex
    for i in instances(Collection)
        @test length(d[i]) == h.npart[Int(i)]
    end

    pos = read_gadget2_pos("pot_acc.format2.gadget2", uGadget2)
    @test length(pos) == 1000
    @test pos[1] == PVector(-0.02657494880259037,
                            -0.040125735104084015,
                            -0.006172948982566595, u"kpc")
end

@testset "GadgetFormat2Fields" begin
    # Units are defined in PhysicalParticles, here we test just to show them
    uAcc = getuAcc(uGadget2)
    uPot = getuEnergyUnit(uGadget2)
    uMass = getuMass(uGadget2)

    @test 1.0*uMass == 1e10u"Msun"
    @test uAcc == u"km^2*kpc^-1*s^-2"
    @test uPot == u"km^2*s^-2"

    h, d = read_gadget2("pot_acc.format2.gadget2", uGadget2, uGadget2)
    @test AstroIO.read_all_mass_from_header(h) == false
    @test AstroIO.read_mass_from_header(h) == [nothing, nothing, nothing, nothing, false, nothing]
    @test d.Acc[1] == PVector(100.51215f0, 146.3331f0, 22.542002f0, uAcc)
    @test d.Potential[1] == -8.796568f0uPot
    @test d.Mass[1] == 1.0f-8*uMass
    @test d.Acc[20] == PVector(1216.8761f0, 868.4943f0, 874.938f0, uAcc)
    @test d.Potential[20] == -37.115158f0uPot
    @test d.Mass[20] == 1.0f-8*uMass

    h, d = read_gadget2("pot_acc.format2.gadget2", uAstro, uGadget2)
    @test d.Acc[20] == PVector(1272.7795f0, 908.3931f0, 915.1328f0, u"kpc*Gyr^-2")
    @test d.Potential[20] == -38.820236f0u"kpc^2*Gyr^-2"
    @test d.Mass[20] == 100.0f0u"Msun"

    write_gadget2_format2("pot_acc.format2.test.gadget2", h, d, acc = true, pot = true)

    h, d = read_gadget2("pot_acc.format2.test.gadget2", uGadget2, uGadget2)
    @test d.Acc[20] == PVector(1216.8761f0,
                               868.4943f0,
                               874.938f0, uAcc)
    @test d.Potential[20] == -37.11515808105469uPot

    @test !iszero(norm(average(d, :Acc)))
    @test !iszero(norm(average(d, :Potential)))
end

@testset "MassFromHeader" begin
    uMass = getuMass(uGadget2)
    h, d = read_gadget2("gadget_no_mass", uAstro)
    @test AstroIO.read_mass_from_header(h) == [true,nothing,nothing,nothing,nothing,nothing]
    @test AstroIO.read_all_mass_from_header(h) == true
    @test d.Mass[1] == 0.10506896f0u"Msun"
    @test Float32(h.mass[1]) == Float32(ustrip(uMass, d.Mass[1]))
end

@testset "FileIO" begin
    h, d = load("gassphere_littleendian.gadget2", uAstro)
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
    header, data = read_gadget2("gassphere_littleendian.gadget2", uAstro)

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

@testset "ConfParser" begin
    ioconfig = loadconfig("conf-jld2.ini")
    data = loadfromconfig(ioconfig)
    @test length(data) == 10

    #TODO Gadget2
    ioconfig = loadconfig("conf-gadget2.ini")
    # header, data = load(ioconfig)
    # @test
end