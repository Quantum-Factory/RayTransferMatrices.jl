import Pkg;
# use a separate Project for building docs; see
# https://discourse.julialang.org/t/psa-use-a-project-for-building-your-docs/14974
Pkg.activate("docs")
Pkg.instantiate()

# to find parent package when executing this file using
# `include("docs/make.jl")`
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

# to find parent package from examples, doctests, etc. (executed in
# working directory docs/build/)
if !("../.." in LOAD_PATH)
    push!(LOAD_PATH, "../..")
end

using Documenter, RayTransferMatrices

makedocs(
    sitename = "RayTransferMatrices",
    format = Documenter.HTML(
        prettyurls = false # put everything in one file
    ),
    # Uncomment the following only when you know what you're
    # doing!  (it will make all doctests pass by overwritting
    # the expected output)
    #doctest = :fix,
    modules = [RayTransferMatrices],
    strict = true
)
# if necessary, repeat the document generation because the example
# code for Symata.jl only works on the second iteration
makedocs(
    sitename = "RayTransferMatrices",
    format = Documenter.HTML(
        prettyurls = false # put everything in one file
    ),
    modules = [RayTransferMatrices],
    strict = true
)

# re-activate parent Project
Pkg.activate(".")
