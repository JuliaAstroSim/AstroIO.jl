# Header
mutable struct HeaderGadget2
    npart::MVector{6,Int32} # gas, halo, disk, bulge, star, blackhole
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
    @MVector zeros(UInt8, 60)
)

function HeaderGadget2(data;
        Counts = count_gadget_types(data),
        time::Float64 = 0.0,
        redshift::Float64 = 0.0,
        counts_total = Counts,
        nfiles::Int64 = 1
    )
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

function read_gadget2_header(filename::AbstractString)
    f = open(filename, "r")

    temp = read(f, Int32)
    if temp == 256
        seekstart(f)
        header = read_gadget2_header(f)
    elseif temp == 8 # Format 2
        name = String(read(f, 4))
        temp1 = read(f, Int32)
        temp2 = read(f, Int32)
        header = read_gadget2_header(f)
    else
        error("Unsupported Gadget2 Format!")
    end

    close(f)

    return header
end

function read_POS!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, uLength, fileuLength)
    temp1 = read(f, Int32)
    Pos = data.Pos
    for i in 1:length(data)
        x = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        y = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        z = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        Pos[i] = PVector(x, y, z)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading positions!\n")
    end
end

function read_VEL!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, uVel, fileuVel)
    temp1 = read(f, Int32)
    Vel = data.Vel
    for i in 1:length(data)
        vx = uconvert(uVel, read(f, Float32) * 1.0 * fileuVel)
        vy = uconvert(uVel, read(f, Float32) * 1.0 * fileuVel)
        vz = uconvert(uVel, read(f, Float32) * 1.0 * fileuVel)
        Vel[i] = PVector(vx, vy, vz)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading velocities!\n")
    end
end

function read_ID!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray)
    temp1 = read(f, Int32)
    ID = data.ID
    for i in 1:length(data)
        ID[i] = Int(read(f, UInt32))
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading IDs!\n")
    end
end

function read_MASS!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, header::HeaderGadget2, uMass, fileuMass)
    read_mass_flag = false
    for i in eachindex(header.npart)
        if header.npart[i] > 0 && header.mass[i] == 0.0
            read_mass_flag = true
            break
        end
    end

    if read_mass_flag
        temp1 = read(f, Int32)
    end

    start = 1
    tail = header.npart[1]
    Mass = data.Mass
    for type in 1:6
        if header.mass[type] == 0.0 # read from file
            for i in start:tail
                Mass[i] = uconvert(uMass, read(f, Float32) * fileuMass)
            end
        else # set using header
            for i in start:tail
                Mass[i] = uconvert(uMass, header.mass[type] * fileuMass)
            end
        end
        start += header.npart[type]
        if type < 6
            tail += header.npart[type+1]
        end
    end
        
    if read_mass_flag
        temp2 = read(f, Int32)
        if temp1 != temp2
            error("Wrong location symbol while reading masses!\n")
        end
    end
end

function read_Entropy!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, NumGas::Int32, uEntropy, fileuEntropy)
    temp1 = read(f, Int32)
    Entropy = data.Entropy
    for i in 1:NumGas
        Entropy[i] = uconvert(uEntropy, read(f, Float32) * 1.0 * fileuEntropy)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Entropy!\n")
    end
end

function read_Density!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, NumGas::Int32, uDensity, fileuDensity)
    temp1 = read(f, Int32)
    Density = data.Density
    for i in 1:NumGas
        Density[i] = uconvert(uDensity, read(f, Float32) * 1.0 * fileuDensity)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Density!\n")
    end
end

function read_HSML!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, NumGas::Int32, uLength, fileuLength)
    temp1 = read(f, Int32)
    HSML = data.Hsml
    for i in 1:NumGas
        HSML[i] = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading Hsml!\n")
    end
end

function read_POT!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, uPot, fileuPot)
    temp1 = read(f, Int32)
    Potential = data.Potential
    for i in 1:length(data)
        Potential[i] = uconvert(uPot, read(f, Float32) * 1.0 * fileuPot)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading potentials\n")
    end
end

function read_ACCE!(f::Union{IOStream,Stream{format"Gadget2"}}, data::StructArray, uAcc, fileuAcc)
    temp1 = read(f, Int32)
    Acc = data.Acc
    for i in 1:length(data)
        accx = uconvert(uAcc, read(f, Float32) * 1.0 * fileuAcc)
        accy = uconvert(uAcc, read(f, Float32) * 1.0 * fileuAcc)
        accz = uconvert(uAcc, read(f, Float32) * 1.0 * fileuAcc)
        Acc[i] = PVector(accx, accy, accz)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading accelerations\n")
    end
end

function read_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units, fileunits;
        acc = false,
        pot = false,
    )
    data = StructArray(Star(units))
    empty!(data)
    for k in 1:6
        append!(data, StructArray(Star(units, collection = GadgetTypes[k]) for i = 1:header.npart[k]))
    end

    if isnothing(units)
        uLength = uVel = uMass = uEntropy = uDensity = uEnergyUnit = uAcc = 1.0
        fileuLength = fileuVel = fileuMass = fileuEntropy = fileuDensity = fileuEnergyUnit = fileuAcc = 1.0
    else
        uLength = getuLength(units)
        uVel = getuVel(units)
        uMass = getuMass(units)
        uEntropy = getuEntropy(units)
        uDensity = getuDensity(units)
        uEnergyUnit = getuEnergyUnit(units)
        uAcc = getuAcc(units)

        fileuLength = getuLength(fileunits)
        fileuVel = getuVel(fileunits)
        fileuMass = getuMass(fileunits)
        fileuEntropy = getuEntropy(fileunits)
        fileuDensity = getuDensity(fileunits)
        fileuEnergyUnit = getuEnergyUnit(fileunits)
        fileuAcc = getuAcc(fileunits)
    end
    
    read_POS!(f, data, uLength, fileuLength)
    read_VEL!(f, data, uVel, fileuVel)
    read_ID!(f, data)
    read_MASS!(f, data, header, uMass, fileuMass)

    # Read Gas Internal Energy Block
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0 && !eof(f)

        if !eof(f)
            read_Entropy!(f, data, NumGas, uEntropy, fileuEntropy)
        end

        if !eof(f)
            read_Density!(f, data, NumGas, uDensity, fileuDensity)
        end

        if !eof(f)
            read_HSML!(f, data, NumGas, uLength, fileuLength)
        end
    end

    if pot
        if eof(f)
            error("No potential data!")
        else
            read_POT!(f, data, uEnergyUnit, fileuEnergyUnit)
        end
    end

    if acc
        if eof(f)
            error("No acceleration data!")
        else
            read_ACCE!(f, data, uAcc, fileuAcc)
        end
    end
    
    return data
end

function read_gadget2_particle_format2(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units, fileunits)
    NumGas = header.npart[1]

    data = StructArray(Star(units))
    empty!(data)
    for k in 1:6
        append!(data, StructArray(Star(units, collection = GadgetTypes[k]) for i = 1:header.npart[k]))
    end

    if isnothing(units)
        uLength = uVel = uMass = uEntropy = uDensity = uEnergyUnit = uAcc = 1.0
        fileuLength = fileuVel = fileuMass = fileuEntropy = fileuDensity = fileuEnergyUnit = fileuAcc = 1.0
    else
        uLength = getuLength(units)
        uVel = getuVel(units)
        uMass = getuMass(units)
        uEntropy = getuEntropy(units)
        uDensity = getuDensity(units)
        uEnergyUnit = getuEnergyUnit(units)
        uAcc = getuAcc(units)

        fileuLength = getuLength(fileunits)
        fileuVel = getuVel(fileunits)
        fileuMass = getuMass(fileunits)
        fileuEntropy = getuEntropy(fileunits)
        fileuDensity = getuDensity(fileunits)
        fileuEnergyUnit = getuEnergyUnit(fileunits)
        fileuAcc = getuAcc(fileunits)
    end

    name = ""
    while !eof(f)
        try
        temp1 = read(f, Int32)
        name = String(read(f, 4))
        temp2 = read(f, Int32)
        skippoint = read(f, Int32)
        catch e
            if isa(e, EOFError)
                return data
            end
        end
        
        if name == "POS "
            read_POS!(f, data, uLength, fileuLength)
        elseif name == "VEL "
            read_VEL!(f, data, uVel, fileuVel)
        elseif name == "ID  "
            read_ID!(f, data)
        elseif name == "MASS"
            read_MASS!(f, data, header, uMass, fileuMass)
        elseif name == "RHO "
            read_Density!(f, data, NumGas, uDensity, fileuDensity)
        elseif name == "HSML"
            read_HSML!(f, data, NumGas, uLength, fileuLength)
        elseif name == "POT "
            read_POT!(f, data, uEnergyUnit, fileuEnergyUnit)
        elseif name == "ACCE"
            read_ACCE!(f, data, uAcc, fileuAcc)
        end
    end

    return data
end

"""
    read_gadget2(filename::AbstractString, units, fileunits = uGadget2; kw...)

Return a Tuple of header and particle data in snapshot file.
`units` is supported by `PhysicalParticles`: `uSI`, `uCGS`, `uAstro`, `uGadget2`, `nothing`.
`fileunits` is the internal units in the file, and will be converted to `units` while reading the file.

# Keywords
- acc::Bool = false : read acceleration data if exist
- pot::Bool = false : read potential data if exist
"""
function read_gadget2(filename::AbstractString, units, fileunits = uGadget2; kw...)
    f = open(filename, "r")

    temp = read(f, Int32)
    if temp == 256
        seekstart(f)
        header = read_gadget2_header(f)
        data = read_gadget2_particle(f, header, units, fileunits; kw...)
    elseif temp == 8 # Format 2
        name = String(read(f, 4))
        temp1 = read(f, Int32)
        temp2 = read(f, Int32)
        header = read_gadget2_header(f)
        data = read_gadget2_particle_format2(f, header, units, fileunits)
    else
        error("Unsupported Gadget2 Format!")
    end

    close(f)

    return header, data
end

function read_gadget2_pos_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units, fileunits)
    NumTotal = sum(header.npart)
    
    if isnothing(units)
        uLength = nothing
        fileuLength = 1.0
    else
        uLength = getuLength(units)
        fileuLength = getuLength(fileunits)
    end
    pos = StructArray(PVector(uLength) for i in 1:NumTotal)

    temp1 = read(f, Int32)
    for i in 1:NumTotal
        x = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        y = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        z = uconvert(uLength, read(f, Float32) * 1.0 * fileuLength)
        pos[i] = PVector(x, y, z)
    end
    temp2 = read(f, Int32)
    if temp1 != temp2
        error("Wrong location symbol while reading positions!\n")
    end

    return pos
end

function read_gadget2_pos_kernel_format2(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units, fileunits)
    while !eof(f)
        temp1 = read(f, Int32)
        name = String(read(f, 4))
        temp2 = read(f, Int32)
        skippoint = read(f, Int32)
        
        if name == "POS "
            return read_gadget2_pos_kernel(f, header, units, fileunits)
        else
            skip(f, temp2)
        end
    end

    @warn "\"POS \" block not found! Returning nothing"
    return nothing
end

"""
    read_gadget2_pos(filename::AbstractString, units, fileunits = uGadget2)

Only read position block and return a StructArray.
"""
function read_gadget2_pos(filename::AbstractString, units, fileunits = uGadget2)
    f = open(filename, "r")

    temp = read(f, Int32)
    if temp == 256
        seekstart(f)
        header = read_gadget2_header(f)
        pos = read_gadget2_pos_kernel(f, header, units, fileunits)
    elseif temp == 8 # Format 2
        name = String(read(f, 4))
        temp1 = read(f, Int32)
        temp2 = read(f, Int32)
        header = read_gadget2_header(f)
        pos = read_gadget2_pos_kernel_format2(f, header, units, fileunits)
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
        Counts[Int(p.Collection)] += 1
    end
    return Counts
end

function count_gadget_types(data::StructArray)
    Counts = MVector{6,Int32}([0,0,0,0,0,0])
    Collection = data.Collection
    for c in Collection
        Counts[Int(c)] += 1
    end
    return Counts
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

function write_POS(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumTotal::Integer, uLength)
    temp = 4 * NumTotal * 3
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(uLength, p.Pos.x)))
                write(f, Float32(ustrip(uLength, p.Pos.y)))
                write(f, Float32(ustrip(uLength, p.Pos.z)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_VEL(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumTotal::Integer, uVel)
    temp = 4 * NumTotal * 3
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(uVel, p.Vel.x)))
                write(f, Float32(ustrip(uVel, p.Vel.y)))
                write(f, Float32(ustrip(uVel, p.Vel.z)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_ID(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumTotal::Integer)
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

function write_MASS_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, header::HeaderGadget2, uMass)
    for i in 1:6
        if header.mass[i] == 0.0
            # if no particle, this would not be executed
            type = GadgetTypes[i]
            for p in data
                if p.Collection == type
                    write(f, Float32(ustrip(uMass, p.Mass)))
                end
            end
        end
    end
end

function write_MASS(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, header::HeaderGadget2, uMass)
    temp = 0
    for i in 1:6
        if header.mass[i] == 0.0 && header.npart[i] != 0
            temp += header.npart[i] * 4
        end
    end

    if temp > 0
        write(f, Int32(temp))
        write_MASS_kernel(f, data, header, uMass)
        write(f, Int32(temp))
    end
end

function write_Entropy(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumGas::Integer, uEntropy)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(uEntropy, p.Entropy)))
        end
    end
    write(f, Int32(temp))
end

function write_Density(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumGas::Integer, uDensity)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(uDensity, p.Density)))
        end
    end
    write(f, Int32(temp))
end

function write_Hsml(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumGas::Integer, uLength)
    temp = NumGas * 4
    write(f, Int32(temp))
    for p in data
        if p.Collection == GAS
            write(f, Float32(ustrip(uLength, p.Hsml)))
        end
    end
    write(f, Int32(temp))
end

function write_POT(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumTotal::Integer, uEnergyUnit)
    temp = 4 * NumTotal
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(uEnergyUnit, p.Potential)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_ACCE(f::Union{IOStream,Stream{format"Gadget2"}}, data::AbstractArray, NumTotal::Integer, uAcc)
    temp = 4 * NumTotal * 3
    write(f, Int32(temp))
    for type in GadgetTypes
        for p in data
            if p.Collection == type
                write(f, Float32(ustrip(uAcc, p.Acc.x)))
                write(f, Float32(ustrip(uAcc, p.Acc.y)))
                write(f, Float32(ustrip(uAcc, p.Acc.z)))
            end
        end
    end
    write(f, Int32(temp))
end

function write_gadget2_particle(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, data, units;
        acc = false,
        pot = false,
    )
    NumTotal = sum(header.npart)

    write_POS(f, data, NumTotal, getuLength(units))
    write_VEL(f, data, NumTotal, getuVel(units))
    write_ID(f, data, NumTotal)
    write_MASS(f, data, header, getuMass(units))

    # Write Gas
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0
        write_Entropy(f, data, NumGas, getuEntropy(units))
        write_Density(f, data, NumGas, getuDensity(units))
        write_Hsml(f, data, NumGas, getuLength(units))
    end

    if pot
        write_POT(f, data, NumTotal, getuEnergyUnit(units))
    end

    if acc
        write_ACCE(f, data, NumTotal, getuAcc(units))
    end

    flush(f)
end

"""
write_gadget2(filename::AbstractString, header::HeaderGadget2, data::AbstractArray, units = uGadget2; format2 = true, kw...)
write_gadget2(filename::AbstractString, data, units = uGadget2; kw...)

Write particle data to file.
`units` are supported by `PhysicalParticles`: `uSI`, `uCGS`, `uAstro`, `uGadget2`, `nothing`.
Gadget2 header is automatedly generated if not provided.

# Keywords
- `format2 = true`: if true, write in format-2 of Gadget2 snapshot
- acc::Bool = false : write acceleration data
- pot::Bool = false : write potential data
"""
function write_gadget2(filename::AbstractString, header::HeaderGadget2, data::AbstractArray, units = uGadget2; format2 = true, kw...)
    if format2
        write_gadget2_format2(filename, header, data, units; kw...)
        return true
    else
        f = open(filename, "w")   
        write_gadget2_header(f, header)
        write_gadget2_particle(f, header, data, units; kw...)
        close(f)
        return true
    end

end

function write_gadget2(filename::AbstractString, data::AbstractArray, units = uGadget2; kw...)
    header = HeaderGadget2(data)
    write_gadget2(filename, header, data, units; kw...)
end

function write_gadget2_format2_block(f::Union{IOStream,Stream{format"Gadget2"}}, name::String, NumBytes::Int64)
    write(f, Int32(8))
    write(f, name)
    write(f, Int32(8 + NumBytes))
    write(f, Int32(8))
end

function write_gadget2_format2(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, data::AbstractArray, units;
        acc = false,
        pot = false,
    )
    write_gadget2_format2_block(f, "HEAD", 256)
    write_gadget2_header(f, header)

    NumTotal = sum(header.npart)

    write_gadget2_format2_block(f, "POS ", 4 * 3 * NumTotal)
    write_POS(f, data, NumTotal, getuLength(units))

    write_gadget2_format2_block(f, "VEL ", 4 * 3 * NumTotal)
    write_VEL(f, data, NumTotal, getuVel(units))

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
        write_MASS_kernel(f, data, header, getuMass(units))
        write(f, Int32(temp))
    end

    # Write Gas
    NumGas = header.npart[1]
    if NumGas > 0 && header.flag_entropy_instead_u > 0
        #TODO Entropy
    
        write_gadget2_format2_block(f, "RHO ", 4 * NumGas)
        write_Density(f, data, NumGas, getuDensity(units))
        
        write_gadget2_format2_block(f, "HSML", 4 * NumGas)
        write_Hsml(f, data, NumGas, getuLength(units))
    end

    if pot
        write_gadget2_format2_block(f, "POT ", 4 * NumTotal)
        write_POT(f, data, NumTotal, getuEnergyUnit(units))
    end

    if acc
        write_gadget2_format2_block(f, "ACCE", 4 * 3 * NumTotal)
        write_ACCE(f, data, NumTotal, getuAcc(units))
    end
end

"""
    write_gadget2_format2(filename::AbstractString, header::HeaderGadget2, data::AbstractArray, units = uGadget2; kw...)
    write_gadget2_format2(filename::AbstractString, data, units = uGadget2; kw...)

Write particle data to format-2 Gadget2 snapshot file.
`units` are supported by `PhysicalParticles`: `uSI`, `uCGS`, `uAstro`, `uGadget2`, `nothing`.
Gadget2 header is automatedly generated if not provided.

# Keywords
- acc::Bool = false : write acceleration data
- pot::Bool = false : write potential data
"""
function write_gadget2_format2(filename::AbstractString, header::HeaderGadget2, data::AbstractArray, units = uGadget2;
        acc = false,
        pot = false,
    )
    f = open(filename, "w")
    write_gadget2_format2(f, header, data, units; acc, pot)
    close(f)
    return true
end

function write_gadget2_format2(filename::AbstractString, data::AbstractArray, units = uGadget2; kw...)
    header = HeaderGadget2(data)
    write_gadget2_format2(filename, header, data, units; kw...)
    return true
end

# FileIO API
import FileIO: Stream, File

function load(s::Stream{format"Gadget2"}, units = uAstro, fileunits = uGadget2)
    header = read_gadget2_header(s)
    data   = read_gadget2_particle(s, header, units, fileunits)
    return header, data
end

function load(f::File{format"Gadget2"}, units = uAstro, fileunits = uGadget2)
    open(f) do s
        header, data = load(s, units, fileunits)
    end
end

function save(s::Stream{format"Gadget2"}, header::HeaderGadget2, data::AbstractArray, units = uGadget2)
    #write_gadget2_header(s, header)
    #write_gadget2_particle(s, header, data)
    write_gadget2_particle(s, header, data, units)
end

function save(f::File{format"Gadget2"}, header::HeaderGadget2, data::AbstractArray, units = uGadget2)
    open(f, "w") do s
        save(s, header, data, units)
    end
end