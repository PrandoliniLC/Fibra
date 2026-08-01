module Fibra

#=
Ho:fiber NALM oscillator — Julia port of the MATLAB script

Units: time [ps], wavelength [nm], freq. [THz], power [W], length [km]
Nonlinear Aeff [µm^2], n2 [1e-16 cm^2/W], gamma [1/(W km)], beta [ps^n/km]

=#
using Printf

include("phys_const.jl")
include("fiber_types.jl")
include("optic_types.jl")
include("math_func.jl")
include("GNLSE.jl")

# =========================================================================
# INPUT FIELD PARAMETERS
# =========================================================================
N2      = 1.0^2                 # Soliton order
tfwhm   = 0.6                   # ps
lamda0  = 2050.0                # pulse central lambda (nm)
f0      = phys_const.c / lamda0            # central pulse frequency (THz)
# -------------------------------------------------------------------------

# =========================================================================
# Numerical parameters
# =========================================================================
nt   = 2^11                      # number of spectral points
time = 70.0                      # ps
dt   = time / nt                 # ps
t    = collect(-time/2 : dt : (time/2 - dt))   # ps

df = 1 / (nt * dt)                             # frequency separation (THz)
f  = collect((-(nt/2)) : 1 : (nt/2 - 1)) .* df # frequency vector (THz)
lambda = phys_const.c ./ (f .+ phys_const.c / lamda0)                # wavelength vector (nm)
w  = 2pi .* f                                  # angular frequency vector

dz  = 0.000001         # longitudinal step (km)
tol = 0.05e-5           # local error tolerance for the propagation solver
# -------------------------------------------------------------------------

# =========================================================================
# Initial fiber parameters
# =========================================================================
smf1 = fiber_types.create_nasmf(fiber_types.PM1950, 1.0, lamda0)

smf2 = smf1
smf2.length = 0.0003

smf3 = smf1
smf3.length = 0.00065

smf4 = smf1
smf4.length = 0.00045

smf5 = smf1
smf5.length = 0.001

amp1 = fiber_types.create_asmf(fiber_types.PM1950, fiber_types.AMFparam, 0.0013, lamda0, lamda0)
# --------------------------------------------------------------------------

# ==========================================================================
# Soliton equations and Initial/Input pulse parameters
# ==========================================================================
println()
P_peak_N = N2 * abs(smf1.beta[3]) / (smf1.gamma * (tfwhm / 1.7627)^2)
E_N      = 2 * abs(smf1.beta[3]) / (smf1.gamma * (tfwhm / 1.7627))
# @printf("N^2=1 soliton: Power = %f [W], Pulse Energy = %f [pJ]\n", P_peak_N, E_N)


P_peak = N2 * abs(smf1.beta[3]) / (smf1.gamma * tfwhm^2) 
u0 = sqrt(P_peak) .* sech.(t ./ tfwhm)          # W^0.5
u0 = math_func.generate_normal_number(0, 10, nt) .* u0
PeakPower = maximum(abs2.(u0))

@show typeof(u0)


@printf("\n----------------------------------------------\n")
@printf("Input Peak Power (W) = %5.2f\n", PeakPower)

b = dt * sum(abs2.(u0))
@printf("Input Pulse Energy in pJ = %5.2f\n", b)
# --------------------------------------------------------------------------

# ==========================================================================
# Filter parameters
# ==========================================================================
filt = optic_types.filter(f0, 30.0, df, 1) # Center Freq., Bandwidth, Freq. Step, order
# --------------------------------------------------------------------------

# ==========================================================================
# coupler parameters
# ==========================================================================
rho     = 0.5    # NALM
rho_out = 0.60   # Output coupler
# --------------------------------------------------------------------------

# ==========================================================================
# NRPS (nonreciprocal phase shift for the NALM)
# ==========================================================================
PhaseShift = 95 * pi / 180
# --------------------------------------------------------------------------

# ==========================================================================
# Being loop
# ==========================================================================
println("\nInteraction Picture Method started")
t_start = time_ns()

spec_z = Vector{Vector{Float64}}()
u_z    = Vector{Vector{ComplexF64}}()

u = copy(u0)
N_trip = 25

#u, _, _ = GNLSE.IP_CQEM_FD(u, dt, dz, smf5, f0, tol, true, true)
u = GNLSE.IP_CQEM_FD(u, dt, dz, amp1, f0, tol, true, true)


#local ufo, ubo, uf, ub, ud, ur, uout   # will hold last-loop values, used after the loop

for ii in 1:N_trip
#    @printf("Round trip %d / %d\n", ii, N_trip)
end
end
