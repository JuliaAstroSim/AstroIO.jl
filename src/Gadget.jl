# Header
mutable struct HeaderGadget2
    npart::MVector{6,Int32} # gas, halo, disk, Bulge, star, blackholw
    mass::MVector{6,Float64}

    time::Float64
    redshift::Float64

    flag_sfr::Int32
    flag_feedback::Int32

    npartTotal::MVector{6,UInt32}

    flag_cooling::Int32

    num_files::Int32

    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64

    flag_stellarage::Int32
    flag_metals::Int32

    npartTotalHighWord::MVector{6,UInt32}

    flag_entropy_instead_u::Int32

    fill_array::MVector{60,UInt8}
end

HeaderGadget2() = HeaderGadget2(MVector{6,Int32}([0,0,0,0,0,0]),
    MVector{6,Float64}([0.0,0.0,0.0,0.0,0.0,0.0]),
    0.0, 0.0, 0, 0,
    MVector{6,UInt32}([0,0,0,0,0,0]),
    0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0,
    MVector{6,UInt32}([0,0,0,0,0,0]), 0,
    @MVector zeros(UInt8, 60))

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

function read_POS!(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, uLength::Units)
    temp1 = read(f, Int32)
    for key in GadgetKeys
        d = data[key]
        for i in 1:length(d)
            x = read(f, Float32) * 1.0 * uLength
            y = read(f, Float32) * 1.0 * uLength
            z = read(f, Float32) * 1.0 * uLength
            d[i] = setproperties!!(d[i], Pos = PVector(x, y, z))
        end
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading positions!\n")
    end
end

function read_VEL!(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, uVel::Units)
    temp1 = read(f, Int32)
    for key in GadgetKeys
        d = data[key]
        for i in 1:length(d)
            vx = uconvert(uVel, read(f, Float32) * 1.0 * u"km/s")
            vy = uconvert(uVel, read(f, Float32) * 1.0 * u"km/s")
            vz = uconvert(uVel, read(f, Float32) * 1.0 * u"km/s")
            d[i] = setproperties!!(d[i], Vel = PVector(vx, vy, vz))
        end
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading velocities!\n")
    end
end

function read_ID!(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict)
    temp1 = read(f, Int32)
    for key in GadgetKeys
        d = data[key]
        for i in 1:length(d)
            d[i] = setproperties!!(d[i], ID = Int(read(f, UInt32)))
        end
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading IDs!\n")
    end
end

function read_MASS!(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, header::HeaderGadget2, uMass::Units)
    read_mass_flag = false
    for i in header.mass
        if i == 0.0
            read_mass_flag = true
            break
        end
    end

    if read_mass_flag
        temp1 = read(f, Int32)
    end

    for type in 1:6
        d = data[GadgetKeys[type]]
        if header.mass[type] == 0.0 # read from file
            for i in 1:length(d)
                d[i] = setproperties!!(d[i], Mass = read(f, Float32) * 1.0e10 * uMass)
            end
        else # set using header
            for i in 1:length(d)
                d[i] = setproperties!!(d[i], Mass = header.mass[type] * 1.0e10 * uMass)
            end
        end
    end
        
    if read_mass_flag
        temp2 = read(f, Int32)
        if temp1 != temp2
            error("Wrong location symbol while reading masses!\n")
        end
    end
end

function read_Entropy!(f::Union{IOStream,Stream{format"Gadget2"}}, d::Array, NumGas::Int32, uEntropy::Units)
    temp1 = read(f, Int32)
    for i in 1:NumGas
        d[i] = setproperties!!(d[i], Entropy = read(f, Float32) * 1.0 * uEntropy)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Entropy!\n")
    end
end

function read_Density!(f::Union{IOStream,Stream{format"Gadget2"}}, d::Array, NumGas::Int32, uDensity::Units)
    temp1 = read(f, Int32)
    for i in 1:NumGas
        d[i] = setproperties!!(d[i], Density = read(f, Float32) * 10e10uDensity)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Density!\n")
    end
end

function read_HSML!(f::Union{IOStream,Stream{format"Gadget2"}}, d::Array, NumGas::Int32, uLength::Units)
    temp1 = read(f, Int32)
    for i in 1:NumGas
        d[i] = setproperties!!(d[i], Hsml = read(f, Float32) * 1.0 * uLength)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Hsml!\n")
    end
end

function read_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units = uAstro)
    data = Dict(
        "gases" => [SPHGas(units, collection = GAS) for i = 1:header.npart[1]],
        "haloes" => [Star(units, collection = HALO) for i = 1:header.npart[2]],
        "disks" => [Star(units, collection = DISK) for i = 1:header.npart[3]],
        "bulges" => [Star(units, collection = BULGE) for i = 1:header.npart[4]],
        "stars" => [Star(units, collection = STAR) for i = 1:header.npart[5]],
        "blackholes" => [Star(units, collection = BLACKHOLE) for i = 1:header.npart[6]],
    )
    
    read_POS!(f, data, getuLength(units))
    read_VEL!(f, data, getuVel(units))
    read_ID!(f, data)
    read_MASS!(f, data, header, getuMass(units))

    # Read Gas Internal Energy Block
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0 && !eof(f)
        d = data["gases"]

        if !eof(f)
            read_Entropy!(f, data["gases"], NumGas, getuEntropy(units))
        end

        if !eof(f)
            read_Density!(f, data["gases"], NumGas, getuDensity(units))
        end

        if !eof(f)
            read_HSML!(f, data["gases"], NumGas, getuLength(units))
        end
    end
    
    return data
end

function read_gadget2_particle_format2(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units = uAstro)
    NumGas = header.npart[1]

    data = Dict(
        "gases" => [SPHGas(units, collection = GAS) for i = 1:header.npart[1]],
        "haloes" => [Star(units, collection = HALO) for i = 1:header.npart[2]],
        "disks" => [Star(units, collection = DISK) for i = 1:header.npart[3]],
        "bulges" => [Star(units, collection = BULGE) for i = 1:header.npart[4]],
        "stars" => [Star(units, collection = STAR) for i = 1:header.npart[5]],
        "blackholes" => [Star(units, collection = BLACKHOLE) for i = 1:header.npart[6]],
    )

    while !eof(f)
        temp1 = read(f, Int32)
        name = String(read(f, 4))
        temp2 = read(f, Int32)
        skippoint = read(f, Int32)
        
        if name == "POS "
            read_POS!(f, data, getuLength(units))
        elseif name == "VEL "
            read_VEL!(f, data, getuVel(units))
        elseif name == "ID  "
            read_ID!(f, data)
        elseif name == "MASS"
            read_MASS!(f, data, header, getuMass(units))
        elseif name == "RHO "
            read_Density!(f, data["gases"], NumGas, getuDensity(units))
        elseif name == "HSML"
            read_HSML!(f, data["gases"], NumGas, getuLength(units))
        end
    end

    return data
end

function read_gadget2(filename::AbstractString, units = uAstro)
    f = open(filename, "r")

    temp = read(f, Int32)
    if temp == 256
        seekstart(f)
        header = read_gadget2_header(f)
        data = read_gadget2_particle(f, header, units)
    elseif temp == 8 # Format 2
        name = String(read(f, 4))
        temp1 = read(f, Int32)
        temp2 = read(f, Int32)
        header = read_gadget2_header(f)
        data = read_gadget2_particle_format2(f, header, units)
    else
        error("Unsupported Gadget2 Format!")
    end

    close(f)

    return header, data
end

function read_gadget2_pos_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units = uAstro)
    NumTotal = sum(header.npart)
    uLength = getuLength(units)
    pos = [PVector(uLength) for i in 1:NumTotal]

    temp1 = read(f, Int32)
    for i in 1:NumTotal
        x = read(f, Float32) * 1.0 * uLength
        y = read(f, Float32) * 1.0 * uLength
        z = read(f, Float32) * 1.0 * uLength
        pos[i] = PVector(x, y, z)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading positions!\n")
    end

    return pos
end

function read_gadget2_pos_kernel_format2(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units = uAstro)
    while !eof(f)
        temp1 = read(f, Int32)
        name = String(read(f, 4))
        temp2 = read(f, Int32)
        skippoint = read(f, Int32)
        
        if name == "POS "
            return read_gadget2_pos_kernel(f, header, units)
        else
            skip(f, temp2)
        end
    end

    @warn "\"POS \" block not found! Returning nothing"
    return nothing
end

function read_gadget2_pos(filename::AbstractString, units = uAstro)
    f = open(filename, "r")

    temp = read(f, Int32)
    if temp == 256
        seekstart(f)
        header = read_gadget2_header(f)
        pos = read_gadget2_pos_kernel(f, header, units)
    elseif temp == 8 # Format 2
        name = String(read(f, 4))
        temp1 = read(f, Int32)
        temp2 = read(f, Int32)
        header = read_gadget2_header(f)
        pos = read_gadget2_pos_kernel_format2(f, header, units)
    else
        error("Unsupported Gadget2 Format!")
    end

    close(f)

    return pos
end

# Write

function count_gadget_types(data::Array{T,N}) where T<:AbstractParticle where N
    Counts = MVector{6,Int32}([0,0,0,0,0,0])
    for p in data
        index = findfirst(x->x==p.Collection, GadgetTypes)
        if !isnothing(index)
            Counts[index] += 1
        end
    end
    return Counts
end

function count_gadget_types(data::Dict)
    Counts = MVector{6,Int32}([0,0,0,0,0,0])
    for type in 1:6
        key = GadgetKeys[type]
        if haskey(data, key)
            Counts[type] = length(data[key])
        end
    end
    return Counts
end

function generate_gadget2_header(data;
                                 Counts = count_gadget_types(data),

                                 time::Float64 = 0.0,
                                 redshift::Float64 = 0.0,
                                 counts_total = Counts,
                                 nfiles::Int64 = 1)
    return HeaderGadget2(
        Counts,
        MVector{6,Float64}([0.0,0.0,0.0,0.0,0.0,0.0]),
        time,
        redshift, 0, 0,
        counts_total,
        0,
        nfiles,

        0.0, 0.3, 0.7, 0.71, 0, 0,
        MVector{6,UInt32}([0,0,0,0,0,0]), 0,
        @MVector zeros(UInt8, 60)
    )
end

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

function write_POS(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumTotal::Integer)
    temp = 4 * NumTotal * 3
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
end

function write_POS(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumTotal::Integer)
    temp = 4 * NumTotal * 3
    write(f, Int32(temp))
    for key in GadgetKeys
        if haskey(data, key)
            for p in data[key]
                write(f, Float32(ustrip(u"kpc", p.Pos.x)))
                write(f, Float32(ustrip(u"kpc", p.Pos.y)))
                write(f, Float32(ustrip(u"kpc", p.Pos.z)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_VEL(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumTotal::Integer)
    temp = 4 * NumTotal * 3
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
end

function write_VEL(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumTotal::Integer)
    temp = 4 * NumTotal * 3
    write(f, Int32(temp))
    for key in GadgetKeys
        if haskey(data, key)
            for p in data[key]
                write(f, Float32(ustrip(u"km/s", p.Vel.x)))
                write(f, Float32(ustrip(u"km/s", p.Vel.y)))
                write(f, Float32(ustrip(u"km/s", p.Vel.z)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_ID(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumTotal::Integer)
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
end

function write_ID(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumTotal::Integer)
    temp = 4 * NumTotal
    write(f, Int32(temp))
    for key in GadgetKeys
        if haskey(data, key)
            for p in data[key]
                write(f, Int32(p.ID))
            end
        end
    end
    write(f, Int32(temp))
end

function write_MASS_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, header::HeaderGadget2)
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
end

function write_MASS_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, header::HeaderGadget2)
    for type in 1:6
        key = GadgetKeys[type]
        if header.mass[type] == 0.0 && haskey(data, key)
            # if no particle, this would not be executed
            for p in data[key]
                write(f, Float32(ustrip(u"Msun", p.Mass) / 1.0e10))
            end
        end
    end
end

function write_MASS(f::Union{IOStream,Stream{format"Gadget2"}}, data, header::HeaderGadget2)
    temp = 0
    for i in 1:6
        if header.mass[i] == 0.0 && header.npart[i] != 0
            temp += header.npart[i] * 4
        end
    end

    if temp > 0
        write(f, Int32(temp))
        write_MASS_kernel(f, data, header)
        write(f, Int32(temp))
    end
end

function write_Entropy(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(u"J/K", p.Entropy)))
        end
    end
    write(f, Int32(temp))
end

function write_Entropy(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    d = data["gases"]
    for p in d
        write(f, Float32(ustrip(u"J/K", p.Entropy)))
    end
    write(f, Int32(temp))
end

function write_Density(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(u"Msun/kpc^3", p.Density) / 1.0e10))
        end
    end
    write(f, Int32(temp))
end

function write_Density(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    d = data["gases"]
    for p in d
        write(f, Float32(ustrip(u"Msun/kpc^3", p.Density) / 1.0e10))
    end
    write(f, Int32(temp))
end

function write_Hsml(f::Union{IOStream,Stream{format"Gadget2"}}, data::Array, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(u"kpc", p.Hsml)))
        end
    end
    write(f, Int32(temp))
end

function write_Hsml(f::Union{IOStream,Stream{format"Gadget2"}}, data::Dict, NumGas::Integer)
    temp = NumGas * 4
    write(f, Int32(temp))
    d = data["gases"]
    for p in d
        write(f, Float32(ustrip(u"kpc", p.Hsml)))
    end
    write(f, Int32(temp))
end

function write_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, data)
    NumTotal = sum(header.npart)

    write_POS(f, data, NumTotal)
    write_VEL(f, data, NumTotal)
    write_ID(f, data, NumTotal)
    write_MASS(f, data, header)

    # Write Gas
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0
        write_Entropy(f, data, NumGas)
        write_Density(f, data, NumGas)
        write_Hsml(f, data, NumGas)
    end
    flush(f)
end

function write_gadget2(filename::AbstractString, header::HeaderGadget2, data)
    f = open(filename, "w")

    write_gadget2_header(f, header)

    write_gadget2_particle(f, header, data)

    close(f)
    return true
end

function write_gadget2(filename::AbstractString, data)
    header = generate_gadget2_header(data)
    write_gadget2(filename, header, data)
end

function write_gadget2_format2_block(f::Union{IOStream,Stream{format"Gadget2"}}, name::String, NumBytes::Int64)
    write(f, Int32(8))
    write(f, name)
    write(f, Int32(8 + NumBytes))
    write(f, Int32(8))
end

function write_gadget2_format2(filename::AbstractString, header::HeaderGadget2, data)
    f = open(filename, "w")

    write_gadget2_format2_block(f, "HEAD", 256)
    write_gadget2_header(f, header)

    NumTotal = sum(header.npart)

    write_gadget2_format2_block(f, "POS ", 4 * 3 * NumTotal)
    write_POS(f, data, NumTotal)

    write_gadget2_format2_block(f, "VEL ", 4 * 3 * NumTotal)
    write_VEL(f, data, NumTotal)

    write_gadget2_format2_block(f, "ID  ", 4 * NumTotal)
    write_ID(f, data, NumTotal)

    # Mass Block
    temp = 0
    for i in 1:6
        if header.mass[i] == 0.0 && header.npart[i] != 0
            temp += header.npart[i] * 4
        end
    end
    if temp > 0
        write_gadget2_format2_block(f, "MASS", temp)
        write(f, Int32(temp))
        write_MASS_kernel(f, data, header)
        write(f, Int32(temp))
    end

    # Write Gas
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0
        #TODO Entropy
    
        write_gadget2_format2_block(f, "RHO ", 4 * NumGas)
        write_Density(f, data, NumGas)
        
        write_gadget2_format2_block(f, "HSML", 4 * NumGas)
        write_Hsml(f, data, NumGas)
    end

    close(f)
    return true
end

function write_gadget2_format2(filename::AbstractString, data)
    header = generate_gadget2_header(data)
    write_gadget2_format2(filename, header, data)
end

# FileIO API
import FileIO: Stream, File

function load(s::Stream{format"Gadget2"}, units = uAstro)
    header = read_gadget2_header(s)
    data   = read_gadget2_particle(s, header, units)
    return header, data
end

function load(f::File{format"Gadget2"})
    open(f) do s
        header, data = load(s)
    end
end

function save(s::Stream{format"Gadget2"}, header::HeaderGadget2, data::Array)
    write_gadget2_header(s, header)
    write_gadget2_particle(s, header, data)
end

function save(f::File{format"Gadget2"}, header::HeaderGadget2, data::Array)
    open(f, "w") do s
        save(s, header, data)
    end
end

function save(f::File{format"Gadget2"}, header::HeaderGadget2, data::Dict)
    open(f, "w") do s
        save(s, header, reduce(vcat, values(data)))
    end
end