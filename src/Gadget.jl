# using PhysicalParticles, Unitful, UnitfulAstro, StaticArrays
# Header
mutable struct HeaderGadget2
    npart::MArray{Tuple{6}, Int32, 1, 6} # gas, halo, disk, Bulge, star, blackholw
    mass::MArray{Tuple{6}, Float64, 1, 6}

    time::Float64
    redshift::Float64

    flag_sfr::Int32
    flag_feedback::Int32

    npartTotal::MArray{Tuple{6}, UInt32, 1, 6}

    flag_cooling::Int32

    num_files::Int32

    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64

    flag_stellarage::Int32
    flag_metals::Int32

    npartTotalHighWord::MArray{Tuple{6}, UInt32, 1, 6}

    flag_entropy_instead_u::Int32

    fill_array::MArray{Tuple{60}, UInt8, 1, 60}
end

HeaderGadget2() = HeaderGadget2(
    MArray{Tuple{6}, Int32}([0,0,0,0,0,0]),
    MArray{Tuple{6}, Float64}([0.0,0.0,0.0,0.0,0.0,0.0]),
    0.0, 0.0, 0, 0,
    MArray{Tuple{6}, UInt32}([0,0,0,0,0,0]),
    0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0,
    MArray{Tuple{6}, UInt32}([0,0,0,0,0,0]), 0,
    @MArray zeros(UInt8, 60)
)

# Read

function read_gadget2_header(f::IOStream)
    Header = HeaderGadget2()

    for i in 1:6
        Header.npart[i] = read(f, Int32)
    end

    for i in 1:6
        Header.mass[i] = read(f, Float64)
    end

    Header.time = read(f, Float64)
    Header.redshift = read(f, Float64)
    Header.flag_sfr = read(f, Int32)
    Header.flag_feedback = read(f, Int32)

    for i in 1:6
        Header.npartTotal[i] = read(f, UInt32)
    end

    Header.flag_cooling = read(f, Int32)
    Header.num_files = read(f, Int32)
    Header.BoxSize = read(f, Float64)
    Header.Omega0 = read(f, Float64)
    Header.OmegaLambda = read(f, Float64)
    Header.HubbleParam = read(f, Float64)
    Header.flag_stellarage = read(f, Int32)
    Header.flag_metals = read(f, Int32)

    for i in 1:6
        Header.npartTotalHighWord[i] = read(f, UInt32)
    end

    Header.flag_entropy_instead_u = read(f, Int32)

    for i = 1:60
        Header.fill_array[i] = read(f, UInt8)
    end

    return Header
end

function read_gadget2(filename::String, showHeader = true)
    f = open(filename, "r")
    @info "Reading data from $filename"

    temp1 = read(f, Int32)
    Header = read_gadget2_header(f)
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol while reading Header!\n")
        quit()
    end
    @info "Header loaded"

    if showHeader
        print(Header)
    end

    NumTotal = sum(Header.npart)

    # Initialize particlles
    gases = [SPHGas() for i=1:Header.npart[1]]
    haloes = [Star() for i=1:Header.npart[2]]
    disks = [Star() for i=1:Header.npart[3]]
    bulges = [Star() for i=1:Header.npart[4]]
    stars = [Star() for i=1:Header.npart[5]]
    blackholes = [Star() for i=1:Header.npart[6]]

    data = Dict(
        :gases => gases,
        :haloes => haloes,
        :disks => disks,
        :bulges => bulges,
        :stars => stars,
        :blackholes => blackholes,
    )
    keys = [:gases, :haloes, :disks, :bulges, :stars, :blackholes]

    # Read positions
    temp1 = read(f, Int32)
    for i in values(data)
        for p in i
            x = read(f, Float32)*u"kpc"
            y = read(f, Float32)*u"kpc"
            z = read(f, Float32)*u"kpc"
            p.Pos = PVector(x, y, z)
        end
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol while reading positions!\n")
        quit()
    end
    @info "Position loaded"

    # Read velocities
    temp1 = read(f, Int32)
    for i in values(data)
        for p in i
            vx = uconvert(u"kpc/Gyr", read(f, Float32)*u"km/s")
            vy = uconvert(u"kpc/Gyr", read(f, Float32)*u"km/s")
            vz = uconvert(u"kpc/Gyr", read(f, Float32)*u"km/s")
            p.Vel = PVector(vx,vy,vz)
        end
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol while reading velocities!\n")
        quit()
    end
    @info "Velocity loaded"

    # Read IDs
    temp1 = read(f, Int32)
    for i in values(data)
        for p in i
            p.ID = read(f, UInt32)
        end
    end
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol while reading IDs!\n")
        quit()
    end
    @info "ID loaded"

    # Read masses
    read_mass_flag = false
    for i in Header.mass
        if i == 0.0
            read_mass_flag = true
            break
        end
    end

    if read_mass_flag
        @info "Reading mass from file"
        temp1 = read(f, Int32)
    end

    for i in 1:6
        @info "test point"
        if Header.mass[i] == 0.0
            for p in data[keys[i]]
                p.Mass = read(f, Float32) * 1.0e10u"Msun"
            end
        else
            for p in data[keys[i]]
                p.Mass = Header.mass[i] * 1.0e10u"Msun"
            end
        end
    end
        
        
    if read_mass_flag
        temp2 = read(f, Int32)
        if temp1!=temp2
            error("Wrong location symbol while reading masses!\n")
            quit()
        end
    end
    @info "Mass loaded"

    # Read Gas Internal Energy Block
    if Header.npart[1] > 0
        # Read Entropy
        temp1 = read(f, Int32)
        for i in values(data)
            for p in i
                p.Entropy = read(f, Float32)*u"J/K"
            end
        end
        temp2 = read(f, Int32)
        if temp1!=temp2
            error("Wrong location symbol while reading Entropy!\n")
            quit()
        end
        @info "Entropy loaded"

        # Read Density
        if !eof(f)
            temp1 = read(f, Int32)
            for i in values(data)
                for p in i
                    p.Density = read(f, Float32)*10e10u"Msun/kpc^3"
                end
            end
            temp2 = read(f, Int32)
            if temp1!=temp2
                error("Wrong location symbol while reading Density!\n")
                quit()
            end
        end
        @info "Density loaded"

        # Read Hsml
        if !eof(f)
            temp1 = read(f, Int32)
            for i in values(data)
                for p in i
                    p.Hsml = read(f, Float32)*u"kpc"
                end
            end
            temp2 = read(f, Int32)
            if temp1!=temp2
                error("Wrong location symbol while reading Hsml!\n")
                quit()
            end
        end
        @info "Hsml loaded"
    end

    close(f)

    return Header, data
end

# Write



# FileIO API