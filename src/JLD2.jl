function write_gadget2_jld(filename::AbstractString, header::HeaderGadget2, data)
    FileIO.save(filename, Dict("header" => header, "data" => data))
    return true
end

function read_gadget2_jld(filename::AbstractString)
    header, data = FileIO.load(filename, "header", "data")
    return header, data
end

function write_jld(filename::AbstractString, data)
    FileIO.save(filename, Dict("data" => data))
    return true
end

function read_jld(filename::AbstractString)
    data = FileIO.load(filename, "data")
    return data
end