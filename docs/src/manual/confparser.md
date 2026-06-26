```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
end
```

# ConfParser

`AstroIO` integrates with [`ConfParser.jl`](https://github.com/jaredsmiller/ConfParser.jl)
so that a single INI file can describe a batch I/O workflow. The dispatch
target is selected by the `format` field, yielding either
[`IOConfigGadget2`](@ref) or [`IOConfigJLD2`](@ref).

## INI file format

```ini
; example.jld2.ini
[io]
filename = ./snapshot.jld2
format   = jld2
label    = data
```

```ini
; example.gadget2.ini
[io]
filename   = ./snapshot.gadget2
format     = gadget2
format2    = true
iofields   = POS, ID, VEL, MASS, RHO, ENTR, HSML, POT, U, TEMP
```

The `iofields` are right-padded with spaces to fit Gadget-2's 4-character
block labels. Up to 4 characters are allowed; anything longer is rejected
with an `AssertionError`.

## Loading

```julia
cfg = loadconfig("example.jld2.ini")
```

returns an `IOConfigJLD2` (or `IOConfigGadget2`). To actually read the data:

```julia
data = loadfromconfig(cfg)
```

`loadfromconfig` for JLD2 simply dispatches to `FileIO.load(filename, label)`;
the Gadget-2 path is currently a placeholder.

## Extending

Add new config structs alongside `IOConfigGadget2` / `IOConfigJLD2` in
`src/ConfParser.jl` and extend the `if/elseif` chain in `loadconfig` to
dispatch on additional `format` values.
