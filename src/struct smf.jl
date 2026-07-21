struct smf
           Aeff::Float64
           n2::Float64
           gamma = 1e4*2*pi*smf1.n2/(lamda0*smf1.Aeff)
           alpha = 0.0
           L::Float64
           betaw                  
end

"""
    MarcatiliMode

Type representing a mode of a hollow capillary as presented in:

Marcatili, E. & Schmeltzer, R.
"Hollow metallic and dielectric waveguides for long distance optical transmission and lasers
(Long distance optical transmission in hollow dielectric and metal circular waveguides,
examining normal mode propagation)."
Bell System Technical Journal 43, 1783–1809 (1964).
"""
struct MarcatiliMode{Ta, Tcore, Tclad, LT} <: AbstractMode
    a::Ta # core radius callable as function of z only, or fixed core radius if a Number
    n::Int # azimuthal mode index
    m::Int # radial mode index
    kind::Symbol # kind of mode (transverse magnetic/electric or hybrid)
    unm::Float64 # mth zero of the nth Bessel function of the first kind
    ϕ::Float64 # overall rotation angle of the mode
    coren::Tcore # callable, returns (possibly complex) core ref index as function of ω
    cladn::Tclad # callable, returns (possibly complex) cladding ref index as function of ω
    model::Symbol # if :full, includes complete influence of complex cladding ref index
    loss::LT # Val{true}() or Val{false}() - whether to include the loss
    aeff_intg::Float64 # Pre-calculated integral fraction for effective area
end

abstract type AbstractField end

"""
    TimeField

Abstract supertype for time-domain only fields.
"""
abstract type TimeField <: AbstractField end

"""
    PulseField(λ0, energy, ϕ, τ0, Itshape)

Represents a temporal pulse with shape defined by `Itshape`.

# Fields
- `λ0::Float64`: the central field wavelength
- `energy::Float64`: the pulse energy
- `power::Float64`: the pulse peak power (**after** applying any spectral phases)
- `ϕ::Vector{Float64}`: spectral phases (CEP, group delay, GDD, TOD, ...)
- `Itshape`: a callable `f(t)` to get the shape of the intensity/power in the time domain
"""
struct PulseField{eT, pT, iT} <: TimeField
    λ0::Float64
    energy::eT
    power::pT
    ϕ::Vector{Float64}
    Itshape::iT
end