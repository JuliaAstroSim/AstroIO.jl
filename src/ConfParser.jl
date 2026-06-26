abstract type IOConfig end

"""
    IOConfigGadget2 <: IOConfig

Configuration object produced by [`loadconfig`](@ref) when the INI file
declares `format = gadget2`. It bundles:

- `filename::String`            — path to the snapshot;
- `format2::Bool`               — whether to emit the Format-2 layout;
- `iofields::Vector{String}`    — list of 4-character Gadget-2 block
  labels (right-padded with spaces) that should be processed.
"""
struct IOConfigGadget2 <: IOConfig
    filename::String
    format2::Bool
    iofields::Vector{String}
end

function Base.show(io::IO, ioconfig::IOConfigGadget2)
    print(io,
        """
        Gadget2 I/O config:
            filename: $(ioconfig.filename)
             format2: $(ioconfig.format2)
            iofields: $(ioconfig.iofields)
        """
    )
end

"""
    loadfromconfig(ioconfig::IOConfigGadget2)

!!! warning "Placeholder"
    Currently returns nothing. Will eventually materialise the particle
    blocks listed in `ioconfig.iofields`.
"""
function loadfromconfig(ioconfig::IOConfigGadget2)
    #TODO how to construct Gadget2Block from iofields

end



"""
    IOConfigJLD2 <: IOConfig

Configuration object produced by [`loadconfig`](@ref) when the INI file
declares `format = jld2`. It bundles:

- `filename::String` — path to the JLD2 file;
- `label::String`    — JLD2 key under which the data is stored.

The companion [`loadfromconfig`](@ref)`(::IOConfigJLD2)` calls
`FileIO.load(filename, label)` to retrieve the data.
"""
struct IOConfigJLD2 <: IOConfig
    filename::String
    label::String
end

function Base.show(io::IO, ioconfig::IOConfigJLD2)
    print(io,
        """
        JLD2 I/O config:
              filename: $(ioconfig.filename)
            data label: $(ioconfig.label)
        """
    )
end

"""
    loadfromconfig(ioconfig::IOConfigJLD2)

Read the JLD2 file referenced by `ioconfig` and return the value stored
under `ioconfig.label`. Thin wrapper around `FileIO.load(filename, label)`.
"""
function loadfromconfig(ioconfig::IOConfigJLD2)
    return FileIO.load(ioconfig.filename, ioconfig.label)
end

"""
    loadconfig(ConfFile::String)

Parse the INI file `ConfFile` and return either an [`IOConfigGadget2`](@ref)
or an [`IOConfigJLD2`](@ref) depending on the `format` field under the
`[io]` section. Use [`loadfromconfig`](@ref) on the returned value to
actually load the data.
"""
function loadconfig(ConfFile::String)
    conf = ConfParse(ConfFile)
    parse_conf!(conf)

    filename = retrieve(conf, "io", "filename")
    format = retrieve(conf, "io", "format")

    if format == "gadget2"
        format2 = retrieve(conf, "io", "format2", Bool)
        iofields = retrieve(conf, "io", "iofields")
        
        # append to 4 Char
        for i in eachindex(iofields)
            Len = length(iofields[i])
            @assert 1 <= Len <= 4 "Field names must have 1~4 characters!"
            iofields[i] *= " "^(4 - Len)
        end

        return IOConfigGadget2(filename, format2, iofields)
    elseif format == "jld2"
        label = retrieve(conf, "io", "label")
        return IOConfigJLD2(filename, label)
    else
        error("Unsupported snapshot format!")
    end
end

