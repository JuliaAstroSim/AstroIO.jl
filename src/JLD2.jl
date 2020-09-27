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

"""
include("AstroIO.jl/src/AstroIO.jl")
using .AstroIO

h,d = read_gadget2("AstroIO.jl/test/gassphere_littleendian.gadget2")

write_jld("temp/jldtest.jld2", d)


data = load("temp/jldtest.jld2", "data")

@load "temp/jldtest.jld2"


write_jld("temp/jldtest.jld2", d)


data = load("temp/jldtest.jld2", "data")

@load "temp/jldtest.jld2"
"""