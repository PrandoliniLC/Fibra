module Fibra

#=
Ho:fiber NALM oscillator — Julia port of the MATLAB script

Units: time [ps], wavelength [nm], freq. [THz], power [W], length [km]
Nonlinear Aeff [µm^2], n2 [1e-16 cm^2/W], gamma [1/(W km)], beta [ps^n/km]

=#
using Printf

include("fiber_types.jl")
include("phys_const.jl")
include("math_func.jl") 

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
# -------------------------------------------------------------------------

# =========================================================================
# Initial fiber parameters
# =========================================================================
smf1 = fiber_types.create_smf(fiber_types.PM1950, 1.0, lamda0)

smf2 = smf1
smf2.length = 0.0003

smf3 = smf1
smf3.length = 0.00065

smf4 = smf1
smf4.length = 0.00045

smf5 = smf1
smf5.length = 0.001

amp1 = fiber_types.create_asmf(fiber_types.PM1950, fiber_types.AMFparam, 0.0013, lamda0)
# --------------------------------------------------------------------------

# ==========================================================================
# Soliton equations and Initial/Input pulse parameters
# ==========================================================================
println()
P_peak_N = N2 * abs(smf1.beta[3]) / (smf1.gamma * (tfwhm / 1.7627)^2)
E_N      = 2 * abs(smf1.beta[3]) / (smf1.gamma * (tfwhm / 1.7627))
@printf("N^2=1 soliton: Power = %f [W], Pulse Energy = %f [pJ]\n", P_peak_N, E_N)


P_peak = N2 * abs(smf1.beta[3]) / (smf1.gamma * tfwhm^2) 
u0 = sqrt(P_peak) .* sech.(t ./ tfwhm)          # W^0.5
u0 = math_func.generate_normal_number(0, 10, nt) .* u0
PeakPower = maximum(abs2.(u0))

@show typeof(u0)


@printf("\n----------------------------------------------\n")
@printf("Input Peak Power (W) = %5.2f\n", PeakPower)

b = dt * sum(abs2.(u0))
@printf("Input Pulse Energy in pJ = %5.2f\n", b)

#=
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


=#

end
