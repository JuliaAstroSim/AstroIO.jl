```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
end
```

# Tools

Two small helpers for renaming files in bulk. Both operate on the current
working directory or an explicit `folder` argument.

## `renamesuffixs`

Replace the extension of every file whose name starts with `namebase` with
`suffix`. Files that already end with `suffix` are left untouched.

```julia
# Rename every file beginning with "test_rename" in the current folder
# to use the .ok extension instead of whatever it has now.
renamesuffixs("./", "test_rename", ".ok")
```

This is useful when collecting output from many simulations that all share
a common prefix.

## `renamereplace`

For every file containing `old` in its name, replace the substring with
`new`. Additional `replace` keyword arguments (e.g. `count = 1`) are
forwarded to `Base.replace`.

```julia
# rename every "before" substring in the folder to "after"
renamereplace("./", "before", "after")
```

## Caveats

- Both functions call `mv(..., force = true)`, so a destination that already
  exists will be overwritten.
- They walk `readdir(folder)` directly and do not recurse. If you need
  recursive behaviour, wrap them in a `walkdir` loop.
