module ABCDBeamTrace

export FreeSpace, Interface, ThinLens, Mirror, Tan, Sag, RTM, Beam
export beamtrace, spotsize, location, discretize
export spotradius, spotradiusfunc
export WithBeam
export transform

using DocStringExtensions
using RecipesBase
import Unitful, Colors, Interpolations

include("elements.jl")
include("beamtrace.jl")
include("comparisons.jl")
include("plot-recipes.jl")

end # module
