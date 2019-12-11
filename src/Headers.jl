mutable struct HeaderGadget2
    npart::SArray{Tuple{6}, Int32, 1, 6} # gas, halo, disk, Bulge, star, blackholw
    mass::SArray{Tuple{6}, Float64, 1, 6}

    time::Float64
    redshift::Float64

    flag_sfr::Int32
    flag_feedback::Int32

    npartTotal::SArray{Tuple{6}, UInt32, 1, 6}

    flag_cooling::Int32

    num_files::Int32

    BoxSize::Float64
    Omega0::Float64
    OmegaLambda::Float64
    HubbleParam::Float64

    flag_stellarage::Int32
    flag_metals::Int32

    npartTotalHighWord::SArray{Tuple{6}, UInt32, 1, 6}

    flag_entropy_instead_u::Int32

    fill_array::SArray{Tuple{60}, UInt8, 1, 60}
end

HeaderGadget2() = HeaderGadget2(SA{Int32}[0,0,0,0,0,0],
                                SA[0.0,0.0,0.0,0.0,0.0,0.0],
                                0.0, 0.0, 0, 0,
                                SA{UInt32}[0,0,0,0,0,0],
                                0, 1, 0.0, 0.3, 0.7, 0.71, 0, 0,
                                SA[0,0,0,0,0,0], 0,
                                @SArray zeros(UInt8, 60))