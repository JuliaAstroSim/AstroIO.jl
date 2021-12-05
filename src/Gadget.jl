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


# abstract type AbstractChemicalComposition end

# mutable struct ChemicalComposition{T<:AbstractFloat}
#     Fe::T
#     Mg::T
# end


struct Gadget2Particle{P, V, A, M, E, F, Et, D, T, dP, dE, Prs, T_1, I<:Integer} <: AbstractParticle3D
    Pos::PVector{P}
    Vel::PVector{V}
    Acc::PVector{A}
    Mass::M
    ID::I
    Collection::Collection

    Ti_endstep::I  # Next integer step on the timeline.
    Ti_begstep::I  # Present integer step on the timeline.
    GravCost::I    # For each two-particle interaction, GravCost += 1

    Potential::E   # Particle potential in the force field
    OldAcc::A      # Save the normalization of acceleration of last step. Useful in Tree n-body method.

    # SPH
    Entropy::Et
    Density::D
    Hsml::P

    Left::F
    Right::F
    NumNgbFound::I

    RotVel::PVector{V}
    DivVel::T_1
    CurlVel::T_1
    dHsmlRho::dP

    Pressure::Prs
    DtEntropy::dE
    MaxSignalVel::V

    Energy::E
    Temperature::T

    # Composition::AbstractChemicalComposition
end

function Gadget2Particle{T, I}(units::Array; id::I = zero(I), collection = STAR) where {T<:AbstractFloat, I<:Integer}
    uLength = getuLength(units)
    uTime = getuTime(units)
    uMass = getuMass(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uEnergy = getuEnergy(units)
    uEnergyUnit = getuEnergyUnit(units)
    uTemperature = getuTemperature(units)
    uEntropy = getuEntropy(units)
    return Gadget2Particle(
        PVector(T, uLength), PVector(T, uVel), PVector(T, uAcc),
        zero(T) * uMass, id, collection,
        zero(I), zero(I), zero(I),
        zero(T) * uEnergyUnit, zero(T) * uAcc,

        zero(T) * uEntropy,
        zero(T) * getuDensity(units), zero(T) * uLength,
        zero(T), zero(T), zero(I),
        PVector(T, uVel),
        zero(T) / uTime, zero(T) / uTime, zero(T) * uLength,
        zero(T) * getuPressure(units),
        zero(T) * uEntropy / uTime,
        zero(T) * uVel,
        zero(T) * uEnergyUnit,
        zero(T) * uTemperature
        # ChemicalComposition(0.0,0.0),
    )
end

function Gadget2Particle{T, I}(::Nothing; id::I = zero(I), collection = STAR) where {T<:AbstractFloat, I<:Integer}
    return Gadget2Particle(
        PVector(T), PVector(T), PVector(T),
        zero(T), id, collection,
        zero(I), zero(I), zero(I),
        zero(T), zero(T),

        zero(T),
        zero(T), zero(T),
        zero(T), zero(T), zero(I),
        PVector(T),
        zero(T), zero(T), zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T)
        # ChemicalComposition(0.0,0.0),
    )
end

# TODO We could consider to change the signature of Gadget2Particle to have only two parametric types.

# Default to Float64, Int64.
function Gadget2Particle(args...; kw...)
    Gadget2Particle{Float64, Int64}(args...; kw... )
end

function Gadget2Particle(T::Type, I::Type, args...; kw...)
    Gadget2Particle{T, I}(args...; kw... )
end


const N_TYPE = length(instances(Collection))

# To be substituted by a config file
const name_mapper = Dict("POS " => :Pos,
                         "ID  " => :ID,
                         "VEL " => :Vel,
                         "ACCE" => :Acc,
                         "MASS" => :Mass,
                         "RHO " => :Density,
                         "ENTR" => :Entropy,
                         "HSML" => :Hsml,
                         "POT " => :Potential,
                         "U   " => :Energy,
                         "TEMP" => :Temperature,
                         )

function get_units(field::Symbol, units::Array)
    units_dict = Dict(:Pos => getuLength,
         :ID => x -> Unitful.NoUnits,
         :Vel => getuVel, 
         :Acc => getuAcc,
         :Mass => getuMass,
         :Density  => getuDensity,
         :Entropy => getuEntropy,
         :Hsml => getuLength,
         :Potential => getuEnergyUnit,
         :Pressure => getuPressure,
         :Temperature => getuTemperature,
         :Energy => getuEnergyUnit,
         )
    return get(units_dict, field, x -> Unitful.NoUnits)(units)
end

# TODO the order here is different from the one in Table 2 of the Gadget2 User Guide
# https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf
const gadget2_format1_names = ["HEAD", "POS ", "VEL ", "ID  ", "MASS",
                               "ENTR", "RHO ", "HSML", "POT ", "ACCE"]

abstract type AbstractGadget2Block end

Gadget2Stream = Union{IOStream,Stream{format"Gadget2"}}

# struct AbstractGadget2Stream{N::Integer} <:Union{IOStream,Stream{format"Gadget2"}}
#     format::N
# end

mutable struct Gadget2Block <: AbstractGadget2Block
    label::Union{String, Missing}
    data_size::Integer
    position::Integer
    partlen::Integer
    data_type::Type
    p_types::AbstractVector{Bool}
end

function Gadget2Block(label::Union{String, Missing}, data_size::Integer, position::Integer)
    return Gadget2Block(label, data_size, position, 0, Number, zeros(Bool, N_TYPE))
end


function Base.show(io::IO, b::Gadget2Block)
    print(io,"$(b.label) - $(b.data_size) ($(b.position))")
end


function set_partlen!(block::Gadget2Block, npart::AbstractVector)
    tot_part = sum(npart)
    # Set the partlen, using our amazing heuristics
    flag = false
    if block.label == "POS " || block.label == "VEL "
        if block.data_size == tot_part * 24
            block.partlen = 24
            block.data_type = Float64
        else
            block.partlen = 12
            block.data_type = Float32
        end
        block.p_types = npart .!= false
        flag = true
    elseif block.label == "ID  "
        # Heuristic for long (64-bit) IDs
        if block.data_size == tot_part * 4
            block.partlen = 4
            block.data_type = Int32
        else
            block.partlen = 8
            block.data_type = Int64
        end
        block.p_types = npart .!= false
        flag = true
    end
    return flag
end

"""
Set up the particle types in the block, with a heuristic,
which assumes that blocks are either fully present or not
for a given particle type
"""
function get_block_one_hot_collection(block::Gadget2Block, npart::AbstractVector)
    tot_part = sum(npart)
    p_types = zeros(Bool, N_TYPE)
    if block.data_size == tot_part * block.partlen
        p_types[:] .= true
        return p_types
    end
    for blocknpart in 1:N_TYPE-1
        # iterate of different possible combinations of particles in the block
        # we stop when we can we match the length of the block
        for perm in permutations([1:N_TYPE...], blocknpart)
            # the 64-bit calculation is important here
            if block.data_size == sum(npart[perm]) * block.partlen
                p_types[perm] .= true
                return p_types
            end
        end
    end
    throw(DomainError(block.label, "Could not determine particle types for block"))
end

"""
Get the dimensionality of the block.
eg, 3 for POS, 1 for most other things
"""
function get_block_dim(block::Gadget2Block)
    return block.partlen ÷ sizeof(block.data_type)
end

function read_mass_from_header(header::HeaderGadget2)
    flag = true
    for i in eachindex(header.npart)
        if header.npart[i] > 0 && header.mass[i] == 0.0
            flag = false
            break
        end
    end
    flag
end

function get_blocks(f::Union{IOStream,Stream{format"Gadget2"}}, format::Int, header::HeaderGadget2)
    blocks = Vector{Gadget2Block}()
    get_block = format == 1 ? get_block_format1 : get_block_format2
    counter = 1
    seekstart(f)
    while !eof(f)
        block = get_block(f)
        if format == 1
            name = gadget2_format1_names[counter]
            counter +=1
            # Do not create the mass block if we're going to read it from header  
            name == "MASS" && read_mass_from_header(header) && continue
            block.label = name
        else
            block.label == "MASS" && read_mass_from_header(header) && continue
        end

        block.label == "HEAD" && continue

        succ = set_partlen!(block, header.npart)
        if !succ
            """
            Figure out what particles are here and what types
            they have. This also is a heuristic, which assumes
            that blocks are either fully present or not for a
            given particle. It also has to try all
            possibilities of dimensions of array and data type.
            """ 
            for (dim, tp) in ((1, Float32), (1, Float64), (3, Float32), (3, Float64))
                try
                    block.data_type = tp
                    block.partlen = sizeof(tp) * dim
                    block.p_types = get_block_one_hot_collection(block, header.npart)
                    succ = true
                    break
                catch e
                    isa(e, DomainError) ? continue : rethrow(e)
                end
            end
        end

        if !succ
            @warn "Encountered a gadget block $(block.label) which could not be interpreted - is it a strange length or data type (length=$(block.data_size))?"
        end

        # get_block_one_hot_collection(block, header.npart)

        push!(blocks, block)
    end
    seekstart(f)
    return blocks
end

function read_mass_block!(f::Gadget2Stream, arr::Array, header::HeaderGadget2, target_units::Units, source_units::Units)
    start = 1
    tail = header.npart[1]
    for type in 1:6
        if header.mass[type] == 0.0 # read from file
            for i in start:tail
                arr[i] = uconvert(target_units, read(f, Float32) * source_units)
            end
        else # set using header
            for i in start:tail
                arr[i] = uconvert(target_units, header.mass[type] * source_units)
            end
        end
        start += header.npart[type]
        if type < 6
            tail += header.npart[type+1]
        end
    end
end


function read_block!(f::Gadget2Stream, b::Gadget2Block, data::StructArray, header::HeaderGadget2, units::Nothing, ::Array)
    read_block!(f, b, data, header, units, nothing)
end

function read_block!(f::Gadget2Stream, b::Gadget2Block, data::StructArray, header::HeaderGadget2, units::Union{Array, Nothing}, fileunits::Union{Array, Nothing})
    name = b.label
    qty = get(name_mapper, name, nothing)
    if isnothing(qty)
        @warn "Cannot map \"$name \" block to any array in data. Skipping"
        return
    end
    arr = getproperty(data, qty)
    target_units = isnothing(units) ? NoUnits : get_units(qty, units)
    source_units = isnothing(fileunits) ? NoUnits : get_units(qty, fileunits)
    T = b.data_type
    seek(f, b.position)
    if name == "MASS"
        read_mass_block!(f, arr, header, target_units, source_units)
        return
    end
    dim = get_block_dim(b)
    if dim == 1
        for i in 1:length(data)
            arr[i] = uconvert(target_units, read(f, T) * source_units)
        end
    elseif dim == 3
        for i in 1:length(data)
            # workaround waiting for https://github.com/JuliaAstroSim/PhysicalParticles.jl/issues/27 to be resolved
            x = uconvert(target_units, read(f, T) * source_units)
            y = uconvert(target_units, read(f, T) * source_units)
            z = uconvert(target_units, read(f, T) * source_units)
            arr[i] = PVector(T(x), T(y), T(z))
        end
    else
        error("Unknown array dimensionality: $dim")
    end
end

function read_gadget2_header(f::Gadget2Stream, format::Integer)
    if format == 1
        return read_gadget2_header(f)
    elseif format == 2
        skip(f, 4*sizeof(Int32))
        return read_gadget2_header(f)
    else
        error("Unknown Gadget2 file format: $format")
    end
end

function get_data_types(blocks::Vector{Gadget2Block})
    float_type = Type{Any}
    integer_type = Type{Any}
    for block in blocks
        if block.label == "POS "
            float_type = block.data_type
        elseif block.label == "ID  "
            integer_type = block.data_type
        end
    end
    types = (float_type=float_type, integer_type=integer_type)
    types
end

"""
    read_gadget2(filename::AbstractString, units, fileunits = uGadget2; kw...)

Return a Tuple of header and particle data in snapshot file.
`units` is supported by `PhysicalParticles`: `uSI`, `uCGS`, `uAstro`, `uGadget2`, `nothing`.
`fileunits` is the internal units in the file, and will be converted to `units` while reading the file.
"""
function read_gadget2(filename::AbstractString, units, fileunits = uGadget2)
    f = open(filename, "r")
    header, data = read_gadget2(f, units, fileunits)
    close(f)
    header, data
end

function read_gadget2(f::Gadget2Stream, units, fileunits = uGadget2)
    format = decide_file_format(f)
    header = read_gadget2_header(f, format)
    blocks = get_blocks(f, format, header)
    dtypes = get_data_types(blocks)
    tot_part = sum(header.npart)
    data = StructArray([AstroIO.Gadget2Particle(dtypes..., units) for i=1:tot_part])
    # Setup collections
    indexes = [0, cumsum(header.npart)...]
    collections = instances(Collection)
    for k in 1:6
        data.Collection[indexes[k]+1:indexes[k+1]] .= collections[k]
    end
    for b in blocks
        read_block!(f, b, data, header, units, fileunits)
    end
    header, data
end

function decide_file_format(f::Gadget2Stream)
    seekstart(f)
    temp = read(f, Int32)
    seekstart(f)
    if temp == 256
        return 1
    elseif temp == 8 # Format 2
        return 2
    else
        error("Unsupported Gadget2 Format!")
    end
end

function get_block_format1(f::Gadget2Stream)
    # Structure of the block:
    # <SIZE><DATA><SIZE>
    block_start = position(f)
    data_size = read(f, Int32)
    block_size = data_size + 2 * sizeof(Int32)
    data_start = block_start + sizeof(Int32)
    seek(f, block_start + data_size + sizeof(Int32))
    data_size_after_block = read(f, Int32)
    @assert data_size == data_size_after_block "data size before/after block data do not match"
    # Move to the end
    seek(f, block_start + block_size)
    return Gadget2Block(missing, data_size, data_start)
end


function get_block_format2(f::Gadget2Stream)
    # Structure of the block:
    # <SKIP><NAME><SIZE+8><SKIP><SIZE><DATA><SIZE>
    block_start = position(f)
    skip(f, sizeof(Int32))
    name = String(read(f, 4))
    block_size = read(f, Int32)
    skip(f, sizeof(Int32))
    data_size = read(f, Int32)
    data_start = position(f)
    @assert data_start == block_start + 20 "Data starting position ($data_start) not 20 bytes away from block start ($block_start)"
    @assert block_size == data_size + 8 "In block $name: size ($block_size) != data_size + 8 ($(data_size+8))"
    # Move to the end
    seek(f, block_start + 20 + data_size)
    data_size_after_block = read(f, Int32)
    @assert data_size == data_size_after_block "In block $name: data size before/after block data do not match"
    return Gadget2Block(name, data_size, data_start)
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

function read_gadget2_pos_kernel(f::Union{IOStream,Stream{format"Gadget2"}}, header::HeaderGadget2, units, fileunits)
    NumTotal = sum(header.npart)
    
    if isnothing(units)
        uLength = 1.0
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
            write(f, Float32(ustrip(uDensity, p.Density) / 1.0e10))
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
    header, data = read_gadget2(s, units, fileunits)
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