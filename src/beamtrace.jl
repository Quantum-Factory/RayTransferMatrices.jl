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
struct Beam{L,CL,N}
    λ::L # vacuum wavelength
    z::L # position along beam axis
    n::N # inde of refraction aka optical density
    x::L # radial extent
    k::N # (angle/sin/tan of) slope of beam
    q::CL # complex beam parameter
end
function Beam(
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
    return Beam{L,CL,N}(
        pλ,
        pz,
        n,
        zero(pλ),
        zero(pn),
        pz - pz0 + 1im * pzR
    )
end
@deprecate Beam(λ_beam, w0_beam) Beam(λ=λ_beam, w0=w0_beam)
@deprecate Beam(λ_beam, w0_beam, n0) Beam(λ=λ_beam, w0=w0_beam, n=n0)

"""

Return `η`, the ratio of optical densities (aka optical indices).

$(SIGNATURES)

The next optical density is in the numerator and the previous optical
density is in the denominator.

"""
η(e::Interface) = e.η
η(e::Union{Tan,Sag}) = η(e.e)
η(e::Element) = 1

"""

Return the effective beam propagation length of an element.

$(SIGNATURES)

The distance is measured along the beam's direction of propagation
(beam axis).

!!! note "To Do"
     The effective propagation length may deviate from the physical
     distance along the beam axis if the optical density deviates from
     unity. Currently, this is not taken into account. Do take it into
     account.

"""
dz(e::FreeSpace) = e.L
dz(e::Union{Tan,Sag}) = dz(e.e)
dz(e::Element{L,N}) where {L,N} = zero(L)

"""

Propagate a [`Beam`](@ref) through an [`Element`](@ref), a System
(represented by a vector of elements), or through a ray transfer
matrix.

$(SIGNATURES)

Note that the `Beam` is the second argument. If a ray transfer matrix
is given as first argument, observe the optional keyword arguments to
pass a propagation distance `dz` and a ratio `η` of new to old
refractive index.

"""
function transform(system::Vector{<:Element}, Γ::Beam)
    for i = 1:length(system)
        Γ = transform(system[1], Γ)
    end
    return Γ
end
transform(e::Element, Γ::Beam) =
    transform(Matrix(e), Γ; dz = dz(e), η = η(e))
transform(m::Matrix, Γ::Beam; dz = zero(m[1,1]), η = one(m[1,1])) =
    Beam(Γ.λ, Γ.z+dz, Γ.n/η, (m*[Γ.x,Γ.k])..., /((m*[Γ.q,1])...))
# the above is rather succinct, hence here an alternative
# implementation for reference
#
# To Do: Remove one of these implementations from the source code
function alt_transform(
    m::Matrix,
    Γ::Beam;
    dz::L = zero(m[1,2]),
    η::N = one(m[1,1])
) where {L,N}
    # new types
    L2 = float(L)
    N2 = float(N)
    CL2 = complex(L2)
    # ray transfer matrix values ABCD
    A::N2 = m[1,1]
    B::L2 = m[1,2]
    C = m[2,1]
    D::N2 = m[2,2]
    # new values
    z::L2 = Γ.z + dz
    xk = m * [Γ.x, Γ.k]
    q2::CL2 = (A*Γ.q + B) / (C*Γ.q + D)
    return Beam{L2,CL2,N2}(
        Γ.λ,
        z,
        (Γ.n / η),
        xk[1],
        xk[2],
        q2
    )
end

"""

Propagate a [`Beam`](@ref) through a system.

$(SIGNATURES)

The first argument, the optical system, is given as a vector of
[`Element`](@ref). The return value is a vector consisting of the
unpropagated beam, and the beams propagated through the first one, the
first two, etc., and, finally, all elements of which the system
consists.

Note that the [`Beam`](@ref) is the second argument.

"""
function beamtrace(elems::Vector{<:Element}, Γ0::Beam)
    Γs = Vector{Beam}(undef, length(elems)+1)
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
spotradius(Γ::Beam) = /(-Γ.λ, π*Γ.n*imag(1/Γ.q)) |> sqrt

@deprecate spotsize(beam) spotradius(beam)

"""

Return the location of a beam, measured along the beam axis.

$(SIGNATURES)

"""
location(Γ::Beam) = Γ.z

"""

Return a function that calculates the beam's waist radius as a
function of the propagation distance.

Expects an optical system (e.g. a vector or collection of
[`Element`](@ref)) as first argument, and a [`Beam`](@ref) as second
argument. The optional keyword argument `outside` can be used to
define a return value for beam positions outside those covered by
the system; the default is to throw a DomainError.

"""
function spotradiusfunc(elements, beam::Beam; outside=nothing)
    beams = Vector{Beam}(undef, length(elements)+1)
    beams[1] = beam
    for (i, el) in enumerate(elements)
        beam = transform(el, beam)
        beams[i+1] = beam
    end
    return function(z)
        if z < beams[1].z || z > beams[length(beams)].z
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
        while beams[i+1].z < z
            i += 1
        end
        i = max(i, 1)
        # add the effect of a FreeSpace element to reach the requested
        # beam position z
        beam = transform(FreeSpace(z - beams[i].z), beams[i])
        return spotradius(beam)
    end
end

"""

Discretize a system by splitting each [`Element`](@ref) that occupies
space.

$(SIGNATURES)

Each element that occupies space is split into `N` appropriately
shortened versions of itself. A vector of elements is returned.

"""
discretize(e::FreeSpace, N::Int) = fill(FreeSpace(e.L/N), N)
discretize(e::Element, N::Int) = e
discretize(els::Vector{<:Element}, N::Int) = vcat(discretize.(els,N)...)
