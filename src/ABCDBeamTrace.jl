module ABCDBeamTrace

export GaussianBeam, GeometricBeam
export FreeSpace, Interface, ThinLens, Mirror, Tan, Sag
export beamtrace, location, discretize
export ior, radialpos, slope
export spotradius, spotradiusfunc
export transform
# deprecated methods
export Beam, spotsize, RTM

using DocStringExtensions
using RecipesBase
import Unitful, Colors, Interpolations

include("elements.jl")
include("beam-abstract.jl")
include("beam-geometric.jl")
include("beam-gaussian.jl")
include("beam-general.jl")
include("comparisons.jl")
include("plot-recipes.jl")

end # module
