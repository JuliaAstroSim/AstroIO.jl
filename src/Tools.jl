function renamesuffixs(folder::AbstractString, namebase::AbstractString, suffix::AbstractString)
    files = readdir(folder)
    for oldname in files
        if startswith(oldname, namebase)
            s = splitext(oldname)
            if last(s) != suffix
                body = first(s)
                mv(joinpath(folder, oldname), joinpath(folder, string(body, suffix)))
            end
        end
    end
end