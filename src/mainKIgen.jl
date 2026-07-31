
using FFTW
using Random
using LinearAlgebra
using Plots
using Printf
gr()

# -----------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------
const c0 = 299792458.0        # speed of light, m/s
const c  = c0 / 1e3           # speed of light in nm/ps

# -----------------------------------------------------------------------
# FiberParams
# -----------------------------------------------------------------------
mutable struct FiberParams
    Aeff::Float64
    n2::Float64
    gamma::Float64
    alpha::Float64
    L::Float64
    betaw::Vector{Float64}
    raman::Int
    ssp::Int
    # gain-medium fields (NaN when unused, mirrors MATLAB's "field not present")
    gssdB::Float64
    PsatdBm::Float64
    lamda_gain::Float64
    landa_bw::Float64
    fc::Float64
    fbw::Float64
end

FiberParams(Aeff, n2, gamma, alpha, L, betaw, raman, ssp) =
    FiberParams(Aeff, n2, gamma, alpha, L, betaw, raman, ssp,
                NaN, NaN, NaN, NaN, NaN, NaN)

mutable struct FilterParams
    lamda_c::Float64
    landa_bw::Float64
    fc::Float64
    f3dB::Float64
    n::Int
end

# -----------------------------------------------------------------------
# Reconstructed helper functions (see porting notes above)
# -----------------------------------------------------------------------

# ASSUMED IMPLEMENTATION — Gaussian random amplitude jitter added to u0
function generate_normal_number(mu::Real, sigma::Real, n::Int)
    return mu .+ sigma .* randn(n)
end

# ASSUMED IMPLEMENTATION — standard lossless 2x2 fiber coupler.
# rho = power splitting ratio (bar port); returns (through, cross)
function coupler(u1, u2, rho::Real)
    t = sqrt(rho)
    k = sqrt(1 - rho)
    uout1 = t .* u1 .+ 1im .* k .* u2
    uout2 = 1im .* k .* u1 .+ t .* u2
    return uout1, uout2
end

# ASSUMED IMPLEMENTATION — super-Gaussian bandpass filter applied in the
# frequency domain. f3dB = 3 dB bandwidth (THz), fc = filter center freq
# (THz), n = super-Gaussian order, f0 = pulse carrier freq, df = freq step.
function filter_gauss(u, f3dB::Real, fc::Real, n::Int, f0::Real, df::Real)
    N = length(u)
    f = ((-N ÷ 2):(N ÷ 2 - 1)) .* df
    H = exp.(-log(sqrt(2)) .* ((f .- (fc - f0)) ./ (f3dB / 2)) .^ (2n))
    U = fftshift(fft(u))
    return ifft(ifftshift(U .* H))
end

# ASSUMED IMPLEMENTATION — full width at half maximum of an intensity
# trace E (real, >=0). Returns (width_in_samples, left_index, right_index).
function fwhm(E::AbstractVector{<:Real})
    Emax, imax = findmax(E)
    half = Emax / 2

    il = imax
    while il > 1 && E[il] > half
        il -= 1
    end
    ir = imax
    while ir < length(E) && E[ir] > half
        ir += 1
    end
    return (ir - il, il, ir)
end

# ASSUMED IMPLEMENTATION — time-domain chirp/phase analysis.
# Returns (delta_w, Eout, width, range) analogous to phase_fit_t.m
# delta_w = instantaneous angular-frequency deviation (chirp), rad/ps
# Eout    = intensity |u|^2
# width   = FWHM in samples
# range   = index range used for the chirp fit (above `thresh` of peak)
function phase_fit_t(u; thresh::Real = 0.05)
    Eout = abs2.(u)
    phase = unwrap(angle.(u))
    delta_w = -[i == 1 || i == length(phase) ? 0.0 :
                (phase[i+1] - phase[i-1]) / 2 for i in eachindex(phase)]
    w, il, ir = fwhm(Eout)
    range = il:ir
    return delta_w, Eout, w, range
end

# ASSUMED IMPLEMENTATION — spectral-domain phase analysis.
# Returns (phase_w, range, Ew) analogous to phase_fit_w.m
function phase_fit_w(u, w, span_factor::Real)
    U = fftshift(fft(u))
    Ew = abs2.(U)
    Ew = Ew ./ maximum(Ew)
    width, il, ir = fwhm(Ew)
    center = (il + ir) ÷ 2
    half = round(Int, span_factor * width / 2)
    lo = max(1, center - half)
    hi = min(length(Ew), center + half)
    range = lo:hi
    phase_w = unwrap(angle.(U))
    return phase_w, range, Ew
end

# Simple phase-unwrap (Julia's DSP.jl has one built in; reimplemented
# here to avoid adding a dependency)
function unwrap(p::AbstractVector{<:Real}; tol::Real = pi)
    q = copy(p)
    for i in 2:length(q)
        d = q[i] - q[i-1]
        while d > tol
            q[i:end] .-= 2pi
            d = q[i] - q[i-1]
        end
        while d < -tol
            q[i:end] .+= 2pi
            d = q[i] - q[i-1]
        end
    end
    return q
end

# -----------------------------------------------------------------------
# STUB — the actual GNLSE propagation solver. NOT reimplemented.
# -----------------------------------------------------------------------
"""
    IP_CQEM_FD(u, dt, dz, fiber, f0, tol, ...)

Interaction-picture method with adaptive step size (conservation-quantity
error method), propagating `u` through `fiber` over its length `fiber.L`.

STUB ONLY — port your original IP_CQEM_FD.m here. This function is the
numerical core of the whole simulation (handles dispersion via `betaw`,
Kerr nonlinearity via `gamma`, Raman response, self-steepening, and,
for `amf1`, saturable gain). Everything else in this file assumes it
returns `(u_out, n_fft_evals, plotdata)` as in the MATLAB version.
"""
function IP_CQEM_FD(u, dt, dz, fiber::FiberParams, f0, tol, args...)
    error("IP_CQEM_FD is not implemented — port this from your original " *
          "MATLAB source (see porting notes at the top of this file).")
end

# =========================================================================
# INPUT FIELD PARAMETERS
# =========================================================================
N2      = 1.0^2                 # Soliton order
tfwhm   = 0.6                   # ps
lamda0  = 2050.0                # pulse central lambda (nm)
f0      = c / lamda0            # central pulse frequency (THz)

# PM-SM polarization-maintaining single-mode fiber (Nufern PM1950)
smf1 = FiberParams(
    54.929,                             # Aeff (µm^2)
    30.0,                                # n2 (1e-16 cm^2/W)
    0.0,                                  # gamma placeholder, set below
    0.0,                                  # alpha (km^-1)
    0.00035,                              # L (km)
    [0.0, 0.0, -95.0, 0.0],                # betaw (ps^n/km)
    0,                                     # raman off
    0,                                     # self-steepening off
)
smf1.gamma = 1e4 * 2pi * smf1.n2 / (lamda0 * smf1.Aeff)   # W^-1 km^-1

# amf1 gain module
amf1              = deepcopy(smf1)
amf1.L            = 0.0013
amf1.gssdB        = 41.0
amf1.PsatdBm      = 20.0
amf1.lamda_gain   = lamda0
amf1.landa_bw     = 20.0
amf1.fc           = c / amf1.lamda_gain
amf1.fbw          = c / (amf1.lamda_gain)^2 * amf1.landa_bw

smf2 = deepcopy(smf1); smf2.L = 0.0003
smf3 = deepcopy(smf1); smf3.L = 0.00065
smf4 = deepcopy(smf1); smf4.L = 0.00045
smf5 = deepcopy(smf1); smf5.L = 0.001
smf6 = deepcopy(smf1); smf6.L = 0.0003

# filter parameters
filt = FilterParams(lamda0, 30.0, 0.0, 0.0, 1)
filt.fc   = c / filt.lamda_c
filt.f3dB = c / (filt.lamda_c)^2 * filt.landa_bw

# coupler parameters
rho     = 0.5    # NALM
rho_out = 0.60   # Output coupler

# NRPS (nonreciprocal phase shift for the NALM)
PhaseShift = 95 * pi / 180

# -------------------------------------------------------------------
# Soliton equations
# -------------------------------------------------------------------
println()
P_peak_N = N2 * abs(smf1.betaw[3]) / (smf1.gamma * (tfwhm / 1.7627)^2)
E_N      = 2 * abs(smf1.betaw[3]) / (smf1.gamma * (tfwhm / 1.7627))
@printf("N^2=1 soliton: Power = %f [W], Pulse Energy = %f [pJ]\n", P_peak_N, E_N)

# -------------------------------------------------------------------
# Numerical parameters
# -------------------------------------------------------------------
nt   = 2^11                      # number of spectral points
time = 70.0                      # ps
dt   = time / nt                 # ps
t    = collect(-time/2 : dt : (time/2 - dt))   # ps

df = 1 / (nt * dt)                             # frequency separation (THz)
f  = collect((-(nt/2)) : 1 : (nt/2 - 1)) .* df # frequency vector (THz)
lambda = c ./ (f .+ c / lamda0)                # wavelength vector (nm)
w  = 2pi .* f                                  # angular frequency vector

dz  = 0.000001         # longitudinal step (km)
tol = 0.05e-5           # local error tolerance for the propagation solver

# -------------------------------------------------------------------
# Input field
# -------------------------------------------------------------------
P_peak = N2 * abs(smf1.betaw[3]) / (smf1.gamma * tfwhm^2) / 20
u0 = sqrt(P_peak) .* sech.(t ./ tfwhm)          # W^0.5
u0 = generate_normal_number(0, 10, nt) .* u0
PeakPower = maximum(abs2.(u0))

@printf("\n----------------------------------------------\n")
@printf("Input Peak Power (W) = %5.2f\n", PeakPower)

b = dt * sum(abs2.(u0))
@printf("Input Pulse Energy in pJ = %5.2f\n", b)

# =========================================================================
# PROPAGATE — numerical solution
# =========================================================================
println("\nInteraction Picture Method started")
t_start = time_ns()

spec_z = Vector{Vector{Float64}}()
u_z    = Vector{Vector{ComplexF64}}()

u = copy(u0)
N_trip = 25

local ufo, ubo, uf, ub, ud, ur, uout   # will hold last-loop values, used after the loop

for ii in 1:N_trip
    @printf("Round trip %d / %d\n", ii, N_trip)

    # mirror - filter - smf6 - coupler - smf5
    uf = filter_gauss(u, filt.f3dB, filt.fc, filt.n, f0, df)
    u = uf

    u, _, _ = IP_CQEM_FD(u, dt, dz, smf6, f0, tol, 1, 1)
    u, ud   = coupler(u, 0, rho_out)
    u, _, _ = IP_CQEM_FD(u, dt, dz, smf5, f0, tol, 1, 1)
    uf, ub  = coupler(u, 0, rho)

    # forward light: smf4 - NRPS - smf3 - smf2 - amf1 - smf1
    ufo, _, _ = IP_CQEM_FD(uf, dt, dz, smf4, f0, tol, 1, 1)
    ufo = ufo .* exp(1im * PhaseShift)
    ufo, _, _ = IP_CQEM_FD(ufo, dt, dz, smf3, f0, tol, 1, 1)
    ufo, _, _ = IP_CQEM_FD(ufo, dt, dz, smf2, f0, tol, 1, 1)
    ufo, _, _ = IP_CQEM_FD(ufo, dt, dz, amf1, f0, tol, 1, 1)
    ufo, _, _ = IP_CQEM_FD(ufo, dt, dz, smf1, f0, tol, 1, 1)

    # backward light: smf1 - amf1 - smf2 - smf3 - NRPS - smf4
    ubo, _, _ = IP_CQEM_FD(ub, dt, dz, smf1, f0, tol, 1, 1)
    ubo, _, _ = IP_CQEM_FD(ubo, dt, dz, amf1, f0, tol, 1, 1)
    ubo, _, _ = IP_CQEM_FD(ubo, dt, dz, smf2, f0, tol, 1, 1)
    ubo, _, _ = IP_CQEM_FD(ubo, dt, dz, smf3, f0, tol, 1, 1)
    ubo = ubo .* exp(-1im * PhaseShift)
    ubo, _, _ = IP_CQEM_FD(ubo, dt, dz, smf4, f0, tol, 1, 1)

    u, ur = coupler(ubo, ufo, rho)

    # smf5 - coupler - smf6 - filter - mirror
    u, _, _ = IP_CQEM_FD(u, dt, dz, smf5, f0, tol, 1, 1)
    u, uout = coupler(u, 0, rho_out)
    u, _, _ = IP_CQEM_FD(u, dt, dz, smf6, f0, tol, 1, 1)

    uf = filter_gauss(u, filt.f3dB, filt.fc, filt.n, f0, df)
    u = uf

    # ---- plot temporal input/output pulses ----
    p1 = plot(t, abs2.(u0), label = "input", color = :blue)
    plot!(p1, t, abs2.(uout), label = "output", color = :red,
          xlabel = "t (ps)", ylabel = "|u(z,t)|^2 (W)",
          title = "Initial (blue) and Final (red) Pulse Shapes")
    display(p1)

    # ---- plot output spectrum ----
    spec = abs2.(fftshift(fft(uout)))
    specnorm = spec ./ lambda .^ 2
    specnorm = specnorm ./ maximum(specnorm)
    p2 = plot(c ./ (f .+ f0), specnorm, color = :red,
              xlabel = "lambda (nm)", ylabel = "Normalized Spectrum (a.u.)",
              title = "Output Spectrum")
    display(p2)

    push!(spec_z, specnorm)
    push!(u_z, uout)
end

elapsed = (time_ns() - t_start) / 1e9
@printf("\nSimulation lasted (s) = %5.2f\n", elapsed)
println()

@printf("\n----------------------------------------------\n")

# -------------------------------------------------------------------
# Pulse power and energy summary
# -------------------------------------------------------------------
PeakPowerFinal = maximum(real.(uout .* conj.(uout)))
@printf("Final Peak Power (W) = %5.2f\n", PeakPowerFinal)
bo = dt * sum(abs2.(uout))
@printf("Final Pulse Energy in pJ = %5.2f\n", bo)

PeakPowerR = maximum(real.(ur .* conj.(ur)))
@printf("Rejected Peak Power (W) = %5.2f\n", PeakPowerR)
br = dt * sum(real.(ur .* conj.(ur)))
@printf("Rejected Pulse Energy in pJ = %5.2f\n", br)

PeakPowerD = maximum(real.(ud .* conj.(ud)))
@printf("Diagnostic Peak Power (W) = %5.2f\n", PeakPowerD)
bd = dt * sum(real.(ud .* conj.(ud)))
@printf("Diagnostic Pulse Energy in pJ = %5.2f\n", bd)

PeakPowerF = maximum(real.(uf .* conj.(uf)))
@printf("Filter Peak Power (W) = %5.2f\n", PeakPowerF)
bf = dt * sum(real.(uf .* conj.(uf)))
@printf("Filter Pulse Energy in pJ = %5.2f\n", bf)

# =========================================================================
# RESULT PLOTS
# =========================================================================
U_z    = reduce(hcat, u_z)'          # N_trip x nt
Spec_z = reduce(hcat, spec_z)'       # N_trip x nt

p3 = heatmap(t, 1:N_trip, abs2.(U_z),
             xlabel = "t (ps)", ylabel = "round trip",
             title = "Output Optical Field Evolution", color = :viridis)
display(p3)

p4 = heatmap(c ./ (f .+ f0), 1:N_trip, Spec_z,
             xlabel = "lambda (nm)", ylabel = "round trip",
             title = "Output Spectrum Evolution", color = :viridis)
display(p4)

# -------------------------------------------------------------------
# Phase / chirp plots for uout, ur, ud, uf
# -------------------------------------------------------------------
function plot_phase_summary(u, label::String)
    delta_w, Eout, width, range = phase_fit_t(u)
    p = plot(t, Eout, ylabel = "|u(z,t)|^2 (W)",
             xlabel = @sprintf("t (ps), %.1fps", width * dt),
             title = label, label = "intensity")
    plot!(twinx(p), t[range], delta_w[range], color = :red,
          ylabel = "Chirp (THz)", label = "chirp")
    display(p)

    println()
    E0 = abs2.(u)
    width0, _, _ = fwhm(E0)
    @printf("%s Pulse Width (FWHM) in samples = %f\n", label, width0)

    phase_w, wrange, Ew = phase_fit_w(u, w, 1.6)
    pw = plot(f[wrange], Ew[wrange], label = "spectrum")
    plot!(twinx(pw), f[wrange], phase_w[wrange], color = :red, label = "phase")
    title!(pw, label)
    display(pw)
end

plot_phase_summary(uout, "Output")
plot_phase_summary(ur, "Rejection")
plot_phase_summary(ud, "Diagnostic")
plot_phase_summary(uf, "Filter")