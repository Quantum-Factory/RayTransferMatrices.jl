"""

Construct a Beam object for modeling the propagation of a Gaussian
(i.e. typical laser) beam.

$(SIGNATURES)

This type is parametrized on three underlying types, `L` (for lengths
such as the wavelength or the location along the beam axis), `CL` (for
complex-valued lengths such as the beam parameter), and `N` (for
dimensionless numbers such as the optical density and a normalization
constant `k`).

"""
struct GaussianBeam{L,CL,N} <: AbstractBeam{L,CL,N}
    b::GeometricBeam{L,CL,N} # contains all the uncommented fields
                             # below
    λ::L # vacuum wavelength
    #z::L # position along beam axis
    #n::N # inde of refraction aka optical density
    #x::L # radial extent
    #k::N # (angle/sin/tan of) slope of beam
    q::CL # complex beam parameter
end
function GaussianBeam(
    ;
    λ::Number, # vacuum wavelength
    w0::Number, # waist radius
    z0::Number = zero(λ), # waist position along beam axis
    n::Number = 1 # ior: index of refraction (aka optical density)
)
    # determine dimensional compatibility of types using the trick
    # that zero(one(x)) == zero(x) if and only if x is dimensionless
    if zero(one(n)) != zero(n)
        throw(DomainError("n must be dimensionless"))
    end
    if zero(one(w0 / λ)) != zero(w0 / λ)
        throw(DomainError("ratios of lengths must be dimensionless"))
    end
    if zero(one(z0 / λ)) != zero(z0 / λ)
        throw(DomainError("ratios of lengths must be dimensionless"))
    end
    # promote arguments, prepending p for "promoted" to variable names
    pλ, pw0, pz0 = promote(float(λ), float(w0), float(z0))
    pn = float(n)
    # determine values
    pz = zero(pλ) # position along beam axis (zero)
    pzR = π * pn * pw0^2 / pλ # Rayleigh length (aka Rayleigh range)
    # determine types
    L = typeof(pλ)
    CL = complex(L)
    N = typeof(one(pλ))
    return GaussianBeam{L,CL,N}(
        GeometricBeam{L,CL,N}(
            pz,
            n,
            zero(pλ),
            zero(pn)
        ),
        pλ,
        pz - pz0 + 1im * pzR
    )
end
@deprecate Beam(λ_beam, w0_beam) GaussianBeam(λ=λ_beam, w0=w0_beam)
@deprecate Beam(λ_beam, w0_beam, n0) GaussianBeam(
    λ=λ_beam,w0=w0_beam, n=n0
)

"""

Propagate a [`GaussianBeam`](@ref) through an [`Element`](@ref), a System
(represented by a vector of elements), or through a ray transfer
matrix.

$(SIGNATURES)

Note that the `GaussianBeam` is the second argument. If a ray transfer matrix
is given as first argument, observe the optional keyword arguments to
pass a propagation distance `dz` and a ratio `η` of new to old
refractive index.

"""
function transform(bywhat::Any, ::AbstractBeam) end
function transform(system::Vector{<:Element}, Γ::AbstractBeam)
    for i = 1:length(system)
        Γ = transform(system[1], Γ)
    end
    return Γ
end
transform(e::Element, Γ::AbstractBeam) =
    transform(Matrix(e), Γ; dz = dz(e), η = η(e))
transform(m::Matrix, Γ::GaussianBeam; dz = zero(m[1,1]), η = one(m[1,1])) =
    GaussianBeam(transform(m, Γ.b; dz = dz, η = η), Γ.λ, /((m*[Γ.q,1])...))

"""

Propagate a [`GaussianBeam`](@ref) through a system.

$(SIGNATURES)

The first argument, the optical system, is given as a vector of
[`Element`](@ref). The return value is a vector consisting of the
unpropagated beam, and the beams propagated through the first one, the
first two, etc., and, finally, all elements of which the system
consists.

Note that the [`GaussianBeam`](@ref) is the second argument.

"""
function beamtrace(elems::Vector{<:Element}, Γ0::GaussianBeam)
    Γs = Vector{GaussianBeam}(undef, length(elems)+1)
    Γs[1] = Γ0
    for (ind, elem) in enumerate(elems)
        Γs[ind+1] = transform(elem, Γs[ind])
    end
    return Γs
end

"""

Return the ``1/e^2`` radius of a beam at its current location.

$(SIGNATURES)

At the radius returned by this function, the intensity drops to
``1/e^2``.

"""
spotradius(Γ::GaussianBeam) = /(-Γ.λ, π*ior(Γ)*imag(1/Γ.q)) |> sqrt
@deprecate spotsize(beam) spotradius(beam)

"""

Return a function that calculates the beam's waist radius as a
function of the propagation distance.

Expects an optical system (e.g. a vector or collection of
[`Element`](@ref)) as first argument, and a [`GaussianBeam`](@ref) as
second argument. The optional keyword argument `outside` can be used
to define a return value for beam positions outside those covered by
the system; the default is to throw a DomainError.

"""
function spotradiusfunc(elements, beam::GaussianBeam; outside=nothing) 
    beams = beamtrace(elements, beam)
    return function(z)
        if (
            z < location(beams[1]) ||
            z > location(beams[length(beams)])
        )    
            if isnothing(outside)
                throw(DomainError(string(
                    "system does not cover ",
                    "the requested beam position"
                )))
            end
            return outside
        end
        # find last element before requested beam position z; note
        # that this could be accelerated by a more intelligent
        # algorithm since the beams[i].z should be sorted (and this
        # algorithm already depends on this expected ordering).
        i = 0
        while location(beams[i+1]) < z
            i += 1
        end
        i = max(i, 1)
        # add the effect of a FreeSpace element to reach the requested
        # beam position z
        beam = transform(FreeSpace(z - location(beams[i])), beams[i])
        return spotradius(beam)
    end
end
