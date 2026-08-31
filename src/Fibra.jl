module Fibra

#=
Ho:fiber NALM oscillator — Julia port of the MATLAB script

Units: time [ps], wavelength [nm], freq. [THz], power [W], length [km]
Nonlinear Aeff [µm^2], n2 [1e-16 cm^2/W], gamma [1/(W km)], beta [ps^n/km]

=#
using Printf
using Plots
using FFTW

include("phys_const.jl")
include("fiber_types.jl")
include("optic_types.jl")
include("math_func.jl")
include("GNLSE.jl")

# =========================================================================
# INPUT FIELD PARAMETERS
# =========================================================================
N2 = 1.0^2                 # Soliton order
tfwhm = 0.6                   # ps
lamda0 = 2050.0                # pulse central lambda (nm)
f0 = phys_const.c / lamda0            # central pulse frequency (THz)
# -------------------------------------------------------------------------

# =========================================================================
# Numerical parameters
# =========================================================================
nt = 2^8                      # number of spectral points
time = 70.0                      # ps
dt = time / nt                 # ps
t = collect((-time/2):dt:(time/2-dt))   # ps

df = 1 / (nt * dt)                             # frequency separation (THz)
f = collect((-(nt/2)):1:(nt/2-1)) .* df # frequency vector (THz)
lambda = phys_const.c ./ (f .+ phys_const.c / lamda0)                # wavelength vector (nm)
w = 2pi .* f                                  # angular frequency vector

dz = 0.000001         # longitudinal step (km)
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
E_N = 2 * abs(smf1.beta[3]) / (smf1.gamma * (tfwhm / 1.7627))
# @printf("N^2=1 soliton: Power = %f [W], Pulse Energy = %f [pJ]\n", P_peak_N, E_N)


P_peak = N2 * abs(smf1.beta[3]) / (smf1.gamma * tfwhm^2)
u0 = sqrt(P_peak) .* sech.(t ./ tfwhm)          # W^0.5
#u0 = math_func.generate_normal_number(0, 10, nt) .* u0
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
rho = 0.5    # NALM
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
u_z = Vector{Vector{ComplexF64}}()

u = 2*copy(u0)
uout = u
ur = u
ud = u
uf = u

N_trip = 5

#u, _, _ = GNLSE.IP_CQEM_FD(u, dt, dz, smf5, f0, tol, true, true)
# = GNLSE.IP_CQEM_FD(u, dt, dz, amp1, f0, tol, true, true)


#local ufo, ubo, uf, ub, ud, ur, uout   # will hold last-loop values, used after the loop

for ii in 1:N_trip
      @printf("Round trip %d / %d\n", ii, N_trip)
      # ---- plot temporal input/output pulses ----
      p1 = plot(t, abs2.(u0), label="input", color=:blue)
      plot!(p1, t, abs2.(uout), label="output", color=:red,
            xlabel="t (ps)", ylabel="|u(z,t)|^2 (W)",
            title="Initial (blue) and Final (red) Pulse Shapes")
      display(p1)

      # ---- plot output spectrum ----
      spec = abs2.(fftshift(fft(uout)))
      specnorm = spec ./ lambda .^ 2
      specnorm = specnorm ./ maximum(specnorm)
      p2 = plot(phys_const.c ./ (f .+ f0), specnorm, color=:red,
            xlabel="lambda (nm)", ylabel="Normalized Spectrum (a.u.)",
            title="Output Spectrum")
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
@printf("Final Peak Power (W)          = %5.2f\n", PeakPowerFinal)
bo = dt * sum(abs2.(uout))
@printf("Final Pulse Energy in pJ      = %5.2f\n", bo)

PeakPowerR = maximum(real.(ur .* conj.(ur)))
@printf("Rejected Peak Power (W)       = %5.2f\n", PeakPowerR)
br = dt * sum(real.(ur .* conj.(ur)))
@printf("Rejected Pulse Energy in pJ   = %5.2f\n", br)

PeakPowerD = maximum(real.(ud .* conj.(ud)))
@printf("Diagnostic Peak Power (W)     = %5.2f\n", PeakPowerD)
bd = dt * sum(real.(ud .* conj.(ud)))
@printf("Diagnostic Pulse Energy in pJ = %5.2f\n", bd)

PeakPowerF = maximum(real.(uf .* conj.(uf)))
@printf("Filter Peak Power (W)         = %5.2f\n", PeakPowerF)
bf = dt * sum(real.(uf .* conj.(uf)))
@printf("Filter Pulse Energy in pJ     = %5.2f\n", bf)

# =========================================================================
# RESULT PLOTS
# =========================================================================
U_z = reduce(hcat, u_z)'          # N_trip x nt
Spec_z = reduce(hcat, spec_z)'       # N_trip x nt

@show(size(U_z))
@show(size(Spec_z))
@show(typeof(U_z))
@show(typeof(Spec_z))
@show(typeof(t))
@show(size(t))


p3 = heatmap(t, 1:N_trip, abs2.(U_z),
      xlabel="t (ps)", ylabel="round trip",
      title="Output Optical Field Evolution", color=:viridis)
display(p3)

p5 = heatmap(reverse(lambda), 1:N_trip, Spec_z,
      xlabel="lambda (nm)", ylabel="round trip",
      title="Output Spectrum Evolution", color=:viridis)
display(p5)

# -------------------------------------------------------------------
# Phase / chirp plots for uout, ur, ud, uf
# -------------------------------------------------------------------
#function plot_phase_summary(u, label::String)
#      delta_w, Eout, width, range = phase_fit_t(u)
#      p = plot(t, Eout, ylabel="|u(z,t)|^2 (W)",
#            xlabel=@sprintf("t (ps), %.1fps", width * dt),
#            title=label, label="intensity")
#      plot!(twinx(p), t[range], delta_w[range], color=:red,
#            ylabel="Chirp (THz)", label="chirp")
#      display(p)###

#      println()
#      E0 = abs2.(u)
#      width0, _, _ = fwhm(E0)
#      @printf("%s Pulse Width (FWHM) in samples = %f\n", label, width0)

#      phase_w, wrange, Ew = phase_fit_w(u, w, 1.6)
#      pw = plot(f[wrange], Ew[wrange], label="spectrum")
#      plot!(twinx(pw), f[wrange], phase_w[wrange], color=:red, label="phase")
#      title!(pw, label)
#      display(pw)
#end

#plot_phase_summary(uout, "Output")
#plot_phase_summary(ur, "Rejection")
#plot_phase_summary(ud, "Diagnostic")
#plot_phase_summary(uf, "Filter")
end
