function renamesuffixs(folder::AbstractString, namebase::AbstractString, suffix::AbstractString)
    files = readdir(folder)
    for oldname in files
        if startswith(oldname, namebase)
            body = first(splitext(oldname))
            mv(joinpath(folder, oldname), joinpath(folder, string(body, suffix)))
        end
    end
end