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

    # read_gadget2_header from filename (format 1)
    h2 = AstroIO.read_gadget2_header("gassphere_littleendian.gadget2")
    @test h2.npart == h.npart
    @test h2.mass == h.mass

    # read_any_mass_from_header
    @test AstroIO.read_any_mass_from_header(h) == false
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

    # read_gadget2_header from filename (format 2)
    h2 = AstroIO.read_gadget2_header("pot_acc.format2.gadget2")
    @test h2.npart == h.npart
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
    @test d.Acc[20] == PVector(1272.7795f0, 908.39307f0, 915.13275f0, u"kpc*Gyr^-2")
    @test d.Potential[20] == -38.820236f0u"kpc^2*Gyr^-2"
    @test d.Mass[20] == 100.0f0u"Msun"

    write_gadget2_format2("pot_acc.format2.test.gadget2", h, d, acc = true, pot = true)

    h, d = read_gadget2("pot_acc.format2.test.gadget2", uGadget2, uGadget2)
    @test d.Acc[20] == PVector(1216.8761f0,
                               868.4943f0,
                               874.9379f0, uAcc)
    @test d.Potential[20] == -37.11515808105469uPot

    @test !iszero(norm(average(d, :Acc)))
    @test !iszero(norm(average(d, :Potential)))
end

@testset "MassFromHeader" begin
    uMass = getuMass(uGadget2)
    h, d = read_gadget2("gadget_no_mass", uAstro)
    @test AstroIO.read_mass_from_header(h) == [true,nothing,nothing,nothing,nothing,nothing]
    @test AstroIO.read_all_mass_from_header(h) == true
    @test AstroIO.read_any_mass_from_header(h) == true
    @test d.Mass[1] == 0.10506896f0u"Msun"
    @test Float32(h.mass[1]) == Float32(ustrip(uMass, d.Mass[1]))
end

@testset "Gadget2ParticleConstructors" begin
    # Gadget2Particle with nothing units
    p = AstroIO.Gadget2Particle(nothing; id = 42, collection = STAR)
    @test p.ID == 42
    @test p.Collection == STAR
    @test p.Pos == PVector(0.0, 0.0, 0.0)

    # Gadget2Particle with units
    p2 = AstroIO.Gadget2Particle(uAstro; id = 1, collection = GAS)
    @test p2.ID == 1
    @test p2.Collection == GAS
    @test p2.Pos == PVector(0.0u"kpc", 0.0u"kpc", 0.0u"kpc")

    # Type-parameterized constructor (with explicit Nothing units)
    p3 = AstroIO.Gadget2Particle(Float32, Int32, nothing)
    @test p3 isa AstroIO.Gadget2Particle
    @test p3.ID == 0
    @test typeof(p3.Mass) == Float32

    # Type-parameterized constructor (with units)
    p4 = AstroIO.Gadget2Particle(Float32, Int32, uAstro)
    @test p4 isa AstroIO.Gadget2Particle
    @test p4.Mass == 0.0f0u"Msun"
    @test p4.Pos.x == 0.0f0u"kpc"
end

@testset "GadgetGetUnits" begin
    # get_units for various fields
    @test AstroIO.get_units(:Pos, uAstro) == u"kpc"
    @test AstroIO.get_units(:Vel, uAstro) == u"kpc*Gyr^-1"
    @test AstroIO.get_units(:Mass, uAstro) == u"Msun"
    @test AstroIO.get_units(:ID, uAstro) == Unitful.NoUnits
    @test AstroIO.get_units(:Acc, uAstro) == u"kpc*Gyr^-2"
    @test AstroIO.get_units(:Potential, uAstro) == u"kpc^2*Gyr^-2"
    @test AstroIO.get_units(:Entropy, uAstro) == u"kpc^2*Gyr^-2*Msun*K^-1"
    @test AstroIO.get_units(:Hsml, uAstro) == u"kpc"
    @test AstroIO.get_units(:Density, uAstro) == u"Msun*kpc^-3"
    @test AstroIO.get_units(:Pressure, uAstro) == u"Msun*kpc^-1*Gyr^-2"
    @test AstroIO.get_units(:Temperature, uAstro) == u"K"
    @test AstroIO.get_units(:Energy, uAstro) == u"kpc^2*Gyr^-2"
    # Unknown field returns NoUnits
    @test AstroIO.get_units(:UnknownField, uAstro) == Unitful.NoUnits
end

@testset "GadgetWriteFormat1" begin
    # Test write_gadget2 with format2=false (format 1)
    h, d = read_gadget2("gassphere_littleendian.gadget2", uAstro)
    @test write_gadget2("testFormat1.gadget2", h, d, uAstro; format2 = false)

    # Read back and verify
    h2, d2 = read_gadget2("testFormat1.gadget2", uAstro)
    @test length(d2) == length(d)
    @test h2.npart == h.npart
end

@testset "GadgetWriteFormat2WithData" begin
    # write_gadget2_format2 with auto-generated header from data
    stars = [Star(uAstro) for i = 1:5]
    @test write_gadget2_format2("testFormat2AutoHeader.gadget2", stars, uAstro)

    h, d = read_gadget2("testFormat2AutoHeader.gadget2", uAstro)
    @test length(d) == 5
end

@testset "GadgetHeaderConstructor" begin
    # HeaderGadget2 from data
    stars = [Star(uAstro) for i = 1:5]
    h = AstroIO.HeaderGadget2(stars)
    @test h.npart[Int(STAR)] == 5
    @test sum(h.npart) == 5

    # HeaderGadget2 default constructor
    h0 = AstroIO.HeaderGadget2()
    @test sum(h0.npart) == 0
end

@testset "GadgetCountTypes" begin
    # count_gadget_types on Array of Stars with different collections
    stars = [Star(uAstro; collection = STAR) for i = 1:3]
    halos = [Star(uAstro; collection = HALO) for i = 1:2]
    data = vcat(stars, halos)
    counts = AstroIO.count_gadget_types(data)
    @test counts[Int(STAR)] == 3
    @test counts[Int(HALO)] == 2

    # count_gadget_types on StructArray
    h, d = read_gadget2("gassphere_littleendian.gadget2", uAstro)
    counts2 = AstroIO.count_gadget_types(d)
    @test counts2[Int(GAS)] == 1472
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

    # Verify write_jld can be read back
    d2 = read_jld("testjldGeneral.jld2")
    @test length(d2) == 10
end

@testset "Tools" begin
    # renamesuffixs
    f = open("test_renamesuffixs.test", "w")
    close(f)
    renamesuffixs(pwd(), "test_rename", ".ok")
    @test isfile("test_renamesuffixs.ok")

    # renamereplace
    f = open("test_renamereplace_old.txt", "w")
    close(f)
    renamereplace(pwd(), "_old", "_new")
    @test isfile("test_renamereplace_new.txt")
    @test !isfile("test_renamereplace_old.txt")
end

@testset "PrettyPrint" begin
    h = AstroIO.HeaderGadget2()
    # Test that show doesn't throw
    buf = IOBuffer()
    show(buf, h)
    s = String(take!(buf))
    @test occursin("Gadget2 Header", s)
    @test occursin("Gas", s)
    @test occursin("Star", s)
    @test occursin("BlackHole", s)

    # Gadget2Block show
    block = AstroIO.Gadget2Block("POS ", 100, 0)
    buf2 = IOBuffer()
    show(buf2, block)
    s2 = String(take!(buf2))
    @test occursin("POS", s2)
end

@testset "ConfParser" begin
    ioconfig = loadconfig("conf-jld2.ini")
    data = loadfromconfig(ioconfig)
    @test length(data) == 10

    # IOConfigJLD2 show
    buf = IOBuffer()
    show(buf, ioconfig)
    s = String(take!(buf))
    @test occursin("JLD2 I/O config", s)

    # Gadget2 config loading
    ioconfig_g2 = loadconfig("conf-gadget2.ini")
    @test ioconfig_g2 isa AstroIO.IOConfigGadget2
    @test ioconfig_g2.filename == "./gassphere_littleendian.gadget2"
    @test ioconfig_g2.format2 == true

    # IOConfigGadget2 show
    buf2 = IOBuffer()
    show(buf2, ioconfig_g2)
    s2 = String(take!(buf2))
    @test occursin("Gadget2 I/O config", s2)
end

@testset "Houdini" begin
    stars = [Star(uAstro) for i = 1:5]

    # write_houdini with 3D particles
    write_houdini("testHoudini3D.hcsv", stars, 0.0, uAstro)
    @test isfile("testHoudini3D.hcsv")

    # write_houdini_append
    AstroIO.write_houdini_append("testHoudini3D.hcsv", stars, 1.0, uAstro)
    @test isfile("testHoudini3D.hcsv")

    # write_houdini with AbstractPoint3D (Ball)
    balls = [Ball(uAstro) for i = 1:3]
    write_houdini("testHoudiniBall.hcsv", balls, 0.0, uAstro)
    @test isfile("testHoudiniBall.hcsv")

    # write_houdini_header with Dict
    data_dict = Dict("stars" => stars, "balls" => balls)
    f = open("testHoudiniHeader.hcsv", "w")
    AstroIO.write_houdini_header(f, data_dict)
    close(f)
    @test isfile("testHoudiniHeader.hcsv")

    # write_houdini_data directly with Vector-of-Vectors (bypasses header)
    data_vov = [stars, balls]
    f = open("testHoudiniDataVov.hcsv", "w")
    AstroIO.write_houdini_header(f, stars)  # reuse AbstractParticle3D header
    AstroIO.write_houdini_data(f, data_vov, 0.5, uAstro; pos_ratio = 1.0, vel_ratio = 1.0)
    close(f)
    @test isfile("testHoudiniDataVov.hcsv")

    # NOTE: write_houdini(filename, Dict, ...) is a known source bug:
    # write_houdini_data uses Iterators.flatten on the Dict which yields
    # the keys (Strings) as elements, then accessing p.ID throws
    # `type String has no field ID`. Filed as a known issue.

    # write_houdini_data with pos_ratio and vel_ratio
    f = open("testHoudiniRatio.hcsv", "w")
    AstroIO.write_houdini_header(f, stars)
    AstroIO.write_houdini_data(f, stars, 0.5, uAstro; pos_ratio = 2.0, vel_ratio = 0.5)
    close(f)
    @test isfile("testHoudiniRatio.hcsv")
end
