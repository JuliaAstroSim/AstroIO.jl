"""
    function renamesuffixs(folder::AbstractString, namebase::AbstractString, suffix::AbstractString)

Replace the end of all files whose name contain `namebase` with `suffix`.

## Example

For files `test_1.middle.txt` and `test_2.middle.txt`, with command
```julia
renamesuffixs("./", "test", ".csv")
```
we will get `test_1.middle.csv` and `test_2.middle.csv`
"""
function renamesuffixs(folder::AbstractString, namebase::AbstractString, suffix::AbstractString)
    files = readdir(folder)
    for oldname in files
        if startswith(oldname, namebase)
            s = splitext(oldname)
            if last(s) != suffix
                body = first(s)
                mv(joinpath(folder, oldname), joinpath(folder, string(body, suffix)), force=true)
            end
        end
    end
end

"""
    function renamereplace(folder::AbstractString, old::AbstractString, new::AbstractString; kw...)

Replace string `old` in all files to `new`. Keywords `kw...` are passed to `replace`
"""
function renamereplace(folder::AbstractString, old::AbstractString, new::AbstractString; kw...)
    files = readdir(folder)
    for oldname in files
        if occursin(old, oldname)
            newname = replace(oldname, old=>new; kw...)
            mv(joinpath(folder, oldname), joinpath(folder, newname), force=true)
        end
    end
end