# Header
mutable struct HeaderGadget2
    npart::MArray{Tuple{6},Int32,1,6} # gas, halo, disk, Bulge, star, blackholw
    mass::MArray{Tuple{6},Float64,1,6}

    time::Float64
    redshift::Float64

    flag_sfr::Int32
    flag_feedback::Int32

    npartTotal::MArray{Tuple{6},UInt32,1,6}

    flag_cooling::Int32

    num_files::Int32

    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64

    flag_stellarage::Int32
    flag_metals::Int32

    npartTotalHighWord::MArray{Tuple{6},UInt32,1,6}

    flag_entropy_instead_u::Int32

    fill_array::MArray{Tuple{60},UInt8,1,60}
end

HeaderGadget2() = HeaderGadget2(MArray{Tuple{6},Int32}([0,0,0,0,0,0]),
    MArray{Tuple{6},Float64}([0.0,0.0,0.0,0.0,0.0,0.0]),
    0.0, 0.0, 0, 0,
    MArray{Tuple{6},UInt32}([0,0,0,0,0,0]),
    0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0,
    MArray{Tuple{6},UInt32}([0,0,0,0,0,0]), 0,
    @MArray zeros(UInt8, 60))

# Read

function read_gadget2_header(f::Union{IOStream,Stream{format"Gadget2"}})
    header = HeaderGadget2()

    temp1 = read(f, Int32)

    for i in 1:6
        header.npart[i] = read(f, Int32)
    end

    for i in 1:6
        header.mass[i] = read(f, Float64)
    end

    header.time = read(f, Float64)
    header.redshift = read(f, Float64)
    header.flag_sfr = read(f, Int32)
    header.flag_feedback = read(f, Int32)

    for i in 1:6
        header.npartTotal[i] = read(f, UInt32)
    end

    header.flag_cooling = read(f, Int32)
    header.num_files = read(f, Int32)
    header.BoxSize = read(f, Float64)
    header.Omega0 = read(f, Float64)
    header.OmegaLambda = read(f, Float64)
    header.HubbleParam = read(f, Float64)
    header.flag_stellarage = read(f, Int32)
    header.flag_metals = read(f, Int32)

    for i in 1:6
        header.npartTotalHighWord[i] = read(f, UInt32)
    end

    header.flag_entropy_instead_u = read(f, Int32)

    for i = 1:60
        header.fill_array[i] = read(f, UInt8)
    end

    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Header!\n")
    end
    
    return header
end

function init_data(header::HeaderGadget2, procs = workers())
    gases = [SPHGas() for i = 1:header.npart[1]]
    haloes = [Star() for i = 1:header.npart[2]]
    disks = [Star() for i = 1:header.npart[3]]
    bulges = [Star() for i = 1:header.npart[4]]
    stars = [Star() for i = 1:header.npart[5]]
    blackholes = [Star() for i = 1:header.npart[6]]

    # Initialize particle types
    for p in haloes
        p.Collection = HALO()
    end

    for p in disks
        p.Collection = DISK()
    end

    for p in bulges
        p.Collection = BULGE()
    end

    for p in blackholes
        p.Collection = BLACKHOLE()
    end

    return [gases; haloes; disks; bulges; stars; blackholes]
end

function read_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2)
    NumTotal = sum(header.npart)

    data = init_data(header, procs)

    # Read positions
    temp1 = read(f, Int32)
    for p in data
        x = read(f, Float32) * u"kpc"
        y = read(f, Float32) * u"kpc"
        z = read(f, Float32) * u"kpc"
        p.Pos = PVector(x, y, z)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading positions!\n")
    end
    @info "Position loaded"

    # Read velocities
    temp1 = read(f, Int32)
    for p in data
        vx = uconvert(u"kpc/Gyr", read(f, Float32) * u"km/s")
        vy = uconvert(u"kpc/Gyr", read(f, Float32) * u"km/s")
        vz = uconvert(u"kpc/Gyr", read(f, Float32) * u"km/s")
        p.Vel = PVector(vx, vy, vz)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading velocities!\n")
    end
    @info "Velocity loaded"

    # Read IDs
    temp1 = read(f, Int32)
    for p in data
        p.ID = read(f, UInt32)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading IDs!\n")
    end
    @info "ID loaded"

    # Read masses
    read_mass_flag = false
    for i in header.mass
        if i == 0.0
            read_mass_flag = true
            break
        end
    end

    if read_mass_flag
        @info "Reading mass from file"
        temp1 = read(f, Int32)
    end

    count_temp = 0
    for type in 1:6
        if header.mass[type] == 0.0
            for i in (count_temp + 1):(count_temp + header.npart[type])
                # if no particle, this would not be executed
                data[i].Mass = read(f, Float32) * 1.0e10u"Msun"
            end
        else
            for i in (count_temp + 1):(count_temp + header.npart[type])
                data[i].Mass = header.mass[type] * 1.0e10u"Msun"
            end
        end
    end
        
        
    if read_mass_flag
        temp2 = read(f, Int32)
        if temp1 != temp2
            error("Wrong location symbol while reading masses!\n")
        end
    end
    @info "Mass loaded"

    # Read Gas Internal Energy Block
    NumGas = header.npart[1]
    if header.npart[1] > 0 && header.flag_entropy_instead_u > 0
        # Read Entropy
        temp1 = read(f, Int32)
        for p in data[1:NumGas]
            p.Entropy = read(f, Float32) * u"J/K"
        end
        temp2 = read(f, Int32)
        if temp1 != temp2
            error("Wrong location symbol while reading Entropy!\n")
        end
        @info "Entropy loaded"

        # Read Density
        if !eof(f)
            temp1 = read(f, Int32)
            for p in data[1:NumGas]
                p.Density = read(f, Float32) * 10e10u"Msun/kpc^3"
            end
            temp2 = read(f, Int32)
            if temp1 != temp2
                error("Wrong location symbol while reading Density!\n")
            end
        end
        @info "Density loaded"

        # Read Hsml
        if !eof(f)
            temp1 = read(f, Int32)
            for p in data[1:NumGas]
                p.Hsml = read(f, Float32) * u"kpc"
            end
            temp2 = read(f, Int32)
            if temp1 != temp2
                error("Wrong location symbol while reading Hsml!\n")
            end
        end
        @info "Hsml loaded"
    end
    
    return data
end

function read_gadget2(filename::String)
    f = open(filename, "r")
    @info "Reading data from $filename"
    
    @info "Reading Header"
    header = read_gadget2_header(f)

    @info "Reading Data"
    data = read_gadget2_particle(f, header)

    close(f)

    return header, data
end

# Write

function write_gadget2_header(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2)
    temp = 256

    write(f, Int32(temp))

    for i in header.npart
        write(f, Int32(i))
    end

    for i in header.mass
        write(f, Float64(i))
    end

    write(f, Float64(header.time))
    write(f, Float64(header.redshift))
    write(f, Int32(header.flag_sfr))
    write(f, Int32(header.flag_feedback))

    for i in header.npartTotal
        write(f, UInt32(i))
    end

    write(f, Int32(header.flag_cooling))
    write(f, Int32(header.num_files))
    write(f, Float64(header.BoxSize))
    write(f, Float64(header.Omega0))
    write(f, Float64(header.OmegaLambda))
    write(f, Float64(header.HubbleParam))
    write(f, Int32(header.flag_stellarage))
    write(f, Int32(header.flag_metals))

    for i in header.npartTotalHighWord
        write(f, UInt32(i))
    end

    write(f, Int32(header.flag_entropy_instead_u))

    for i = 1:60
        write(f, UInt8(0))
    end

    write(f, Int32(temp))
    flush(f)
end

function write_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, data::Array)
    NumTotal = sum(header.npart)
    temp = 4 * NumTotal * 3

    @info "  Writing Position"
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(u"kpc", p.Pos.x)))
                write(f, Float32(ustrip(u"kpc", p.Pos.y)))
                write(f, Float32(ustrip(u"kpc", p.Pos.z)))
            end
        end
    end
    write(f, Int32(temp))

    @info "  Writing Velocity"
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(u"km/s", p.Vel.x)))
                write(f, Float32(ustrip(u"km/s", p.Vel.y)))
                write(f, Float32(ustrip(u"km/s", p.Vel.z)))
            end
        end
    end
    write(f, Int32(temp))

    @info "  Writing ID"
    temp = 4 * NumTotal
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Int32(p.ID))
            end
        end
    end
    write(f, Int32(temp))

    @info "  Writing Mass"
    temp = 0
    for i in 1:6
        if header.mass[i] == 0.0 && header.npart[i] != 0
            temp += header.npart[i] * 4
        end
    end

    if temp > 0
        write(f, Int32(temp))

        for i in 1:6
            if header.mass[i] == 0.0
                # if no particle, this would not be executed
                for p in data
                    if p.Collection == GadgetTypes[i]
                        write(f, Float32(ustrip(u"Msun", p.Mass) / 1.0e10))
                    end
                end
            end
        end

        write(f, Int32(temp))
    end

    # Write Gas
    if header.npart[1] > 0 && header.flag_entropy_instead_u > 0
        temp = header.npart[1] * 4
        @info "  Writing Entropy"
        write(f, Int32(temp))
        for p in data
            if p.Collection == GAS()
                write(f, Float32(ustrip(u"J/K", p.Entropy)))
            end
        end
        write(f, Int32(temp))

        @info "  Writing Density"
        write(f, Int32(temp))
        for p in data
            if p.Collection == GAS()
                write(f, Float32(ustrip(u"Msun/kpc^3", p.Density) / 1.0e10))
            end
        end
        write(f, Int32(temp))

        @info "  Writing Hsml"
        write(f, Int32(temp))
        for p in data
            if p.Collection == GAS()
                write(f, Float32(ustrip(u"kpc", p.Hsml)))
            end
        end
        write(f, Int32(temp))
    end
    flush(f)
end

function write_gadget2(filename::String, header::HeaderGadget2, data::Array)
    f = open(filename, "w")

    @info "Writing Header"
    write_gadget2_header(f, header)

    @info "Writing Particle Data"
    write_gadget2_particle(f, header, data)

    close(f)
    @info "Data saved"
    return true
end

# FileIO API

function load(s::Stream{format"Gadget2"})
    header = read_gadget2_header(s)
    data   = read_gadget2_particle(s, header)
    return header, data
end

function load(f::File{format"Gadget2"})
    @info "Loading with FileIO"
    open(f) do s
        header, data = load(s)
    end
end

function save(s::Stream{format"Gadget2"}, header::HeaderGadget2, data::Array)
    write_gadget2_header(s, header)
    write_gadget2_particle(s, header, data)
end

function save(f::File{format"Gadget2"}, header::HeaderGadget2, data::Array)
    @info "Writing with FileIO"
    open(f, "w") do s
        save(s, header, data)
    end
    return true
end