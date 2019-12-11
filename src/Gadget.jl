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

function read_gadget2(filename::String)
    f = open(filename, "r")

    temp1 = read(f, Int32)
    Header = read_gadget2_header(f)
    temp2 = read(f, Int32)
    if temp1!=temp2
        error("Wrong location symbol while reading Header!\n")
        quit()
    end

    NumTotal = sum(Header.npart)

    close(f)

    return Header
end

# Write



# FileIO API