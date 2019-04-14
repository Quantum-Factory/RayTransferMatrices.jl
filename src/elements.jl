"""

An element is the abstract supertype of all elements that together
form an optical system, consistently modeled as a `Vector{<:Element}`.

$(SIGNATURES)

It is parametrized on types L and N, used for numbers of the
respective dimenions length (L) and dimensionless (N, "number").

"""
abstract type Element{L,N} end

"""

An optical element representing propagation over free space.

$(SIGNATURES)

The optical density is assumed to be unity. The propagation length is
the only field in this structure and hence the only argument to the
inherent constructor.

"""
struct FreeSpace{T,N} <: Element{T,N}
    L::T
end
FreeSpace(L) = FreeSpace{typeof(L), typeof(one(L))}(L)

"""

An optical Interface.

$(SIGNATURES)

The Interface has the parameters `η = n1 / n2` where `n1` and `n2` are
the optical densities (aka optical indices) of the previous resp. new
medium, the angle of incidence `aoi`, and a radius of curvature `roc`
which can be infinite. All arguments are optional but the default
arguments have no effect on a beam as the optical densities both
default to zero.

"""
struct Interface{L,N} <: Element{L,N}
    η::N; θ::N; R::L # R > 0 when light hits concave side
end
Interface(; n1=1.0, n2=1.0, η=n2/n1, aoi=0, roc=Inf) =
    Interface(η, aoi, roc)
@deprecate Interface(η,θ) Interface(n1=1,n2=η,aoi=θ)

"""

A thin lens.

$(SIGNATURES)

The parameters are focal length `f` and an optical angle of incidence
`aoi`.

"""
struct ThinLens{L,N} <: Element{L,N}
    f::L; θ::N
end
ThinLens(; f::L, aoi::N = zero(one(f))) where {L,N} =
    ThinLens{L,N}(f, aoi)
@deprecate ThinLens(focallength) ThinLens(f=focallength)

"""

A mirror.

$(SIGNATURES)

The optional keyword arguments are radius of curvature `roc` and angle
of incidence `aoi`.

"""
function Mirror(; roc=Inf, f=0.5*roc, aoi=0)
    if roc ≉ 2f
        throw(ArgumentError("roc and f are incompatible"))
    end
    return ThinLens(f=f, aoi=aoi)
end
@deprecate Mirror(R, θ) Mirror(roc=R, aoi=θ)

"""

Construct a pseudo-element for modeling beam propagation in the
tangential (aka parallel) plane.

$(SIGNATURES)

"""
struct Tan{E,T,N} <: Element{T,N} where {E<:Element{T,N}}
    e::E
end
Tan(e::FreeSpace) = e
Tan(e::Element{T,N}) where {T,N} = Tan{typeof(e),T,N}(e)
Tan(elements::Vector{<:Element}) = Tan.(elements)

"""

Construct a pseudo-element for modeling beam propagation in the
sagittal plane.

$(SIGNATURES)

"""
struct Sag{E,T,N} <: Element{T,N} where {E<:Element{T,N}}
    e::E
end
Sag(e::FreeSpace) = e
Sag(e::Element{T,N}) where {T,N} = Sag{typeof(e),T,N}(e)
Sag(elements::Vector{<:Element}) = Sag.(elements)

"""

`RTM` returns the Ray Transfer (ABCD) matrix associated with the
given, optical [`Element`](@ref) or an entire system (a vector of
`Element`).

$(SIGNATURES)

The matrix is represented as a julia `Matrix` with eltype `number`
unless a specialization applies when using plain numbers (without
units).

"""
RTM(e::FreeSpace) = [1 e.L ; zero(inv(e.L)) 1]
RTM(e::Interface) = [1 zero(e.R) ; (e.η-1)/e.R e.η]
RTM(e::ThinLens) = [1 zero(e.f) ; -inv(e.f) 1]
RTM(e::Tan{ThinLens{T,N},T,N}) where {T,N} =
    RTM(ThinLens(; f = e.e.f * cos(e.e.θ), aoi = e.e.θ))
RTM(e::Sag{ThinLens{T,N},T,N}) where {T,N}  =
    RTM(ThinLens(; f = e.e.f / cos(e.e.θ), aoi = e.e.θ))
# See doi:10.1364/AO.26.000427 for the following matrices
function RTM(e::Tan{Interface{T,N},T,N}) where {T,N}
    θ1, η, R = e.e.θ, e.e.η, e.e.R; θ2 = asin(η*sin(θ1))
    return [cos(θ2)/cos(θ1) zero(L) ;
        (cos(θ2)-η*cos(θ1))/(R*cos(θ1)*cos(θ2)) η*cos(θ1)/cos(θ2)]
end
function RTM(e::Sag{Interface{T,N},T,N}) where {T,N}
    θ1, η, R = e.e.θ, e.e.η, e.e.R; θ2 = asin(η*sin(θ1))
    return [one(N) zero(L) ; (cos(θ2)-η*cos(θ1))/R η]
end
RTM(elements::Vector{<:Element}) = prod(RTM.(elements))
