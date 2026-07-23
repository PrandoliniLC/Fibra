module math_func
using Random


function generate_normal_number(mu, sig, nr)

    Random.seed!(1)

    nor = ones(Float64, nr)

    for i in 1:nr
        u1 = rand()
        u2 = rand()

        z0 = sqrt(-2 * log(u1)) * cos(2π * u2)
        nor[i] = sig * z0 + mu
    end

    return nor
end



"""
    fwhm(x)

Compute the full width at half maximum (FWHM) of a pulse.

# Arguments
- `x`: Vector containing the pulse intensity or amplitude.

# Returns
- `width`: FWHM in samples.
- `I_l`: Left half-maximum index.
- `I_r`: Right half-maximum index.
"""
function fwhm(x::AbstractVector)

    peak, ind_peak = findmax(x)
    half_peak = peak / 2

    # Left and right portions of the pulse
    x_l = @view x[1:ind_peak]
    x_r = @view x[ind_peak+1:end]

    # Distance from the peak to the half-maximum points
    I_l = findfirst(reverse(x_l) .<= half_peak)
    I_r = findfirst(x_r .<= half_peak)

    # Handle the case where no half-maximum crossing is found
    if isnothing(I_l) || isnothing(I_r)
        return nothing, nothing, nothing
    end

    width = I_l + I_r

    # Convert to indices in the original vector
    I_l = ind_peak - I_l
    I_r = ind_peak + I_r

    return width, I_l, I_r
end

end

using DSP          # for unwrap()
using FFTW

"""
    phase_fit_t(u)

Calculate the temporal phase derivative and determine the fitting range
around the pulse peak.

# Arguments
- `u::AbstractVector{<:Complex}` : Complex field.

# Returns
- `delta_w` : First difference of the unwrapped phase.
- `Eu`      : Pulse intensity |u|².
- `width`   : Full Width at Half Maximum (FWHM) in samples.
- `range`   : Indices spanning ±1 FWHM about the pulse peak.
"""
function phase_fit_t(u::AbstractVector{<:Complex})

    # Unwrapped phase
    phase_out = unwrap(angle.(u))

    # Phase derivative
    delta_w = diff(phase_out)

    # Intensity
    Eu = abs2.(u)

    # Peak intensity and location
    uout_max, I_uout_max = findmax(Eu)

    # Full Width at Half Maximum
    width, I_l, I_r = fwhm(Eu)

    # Region around the pulse peak
    first_idx = max(1, floor(Int, I_uout_max - width))
    last_idx  = min(length(Eu), floor(Int, I_uout_max + width))

    range = first_idx:last_idx

    return delta_w, Eu, width, range

end

using FFTW
using DSP               # unwrap()
using Polynomials

"""
    phase_fit_w(u, w, width_fac)

Calculate the spectral phase, perform 3rd- and 4th-order polynomial fits,
and estimate the dispersion coefficients.

# Arguments
- `u::AbstractVector{<:Complex}` : Complex electric field.
- `w::AbstractVector`            : Frequency (or angular frequency) axis.
- `width_fac::Real`              : Width multiplier for the fitting window.

# Returns
- `phase_out` : Unwrapped spectral phase.
- `range`     : Indices used for the polynomial fit.
- `Ew`        : Spectral intensity.
"""
function phase_fit_w(
    u::AbstractVector{<:Complex},
    w::AbstractVector,
    width_fac::Real
)

    # Spectrum
    spec_0 = fftshift(fft(fftshift(u)))

    # Spectral phase
    phase_out = unwrap(angle.(spec_0))

    # Phase derivative (not used later, but retained from MATLAB)
    delta_w = diff(phase_out)

    # Spectral intensity
    Ew = abs2.(spec_0)

    # Peak location
    uout_max, I_uout_max = findmax(Ew)

    # FWHM
    width, I_l, I_r = fwhm(Ew)

    # Fit range
    first_idx = max(1, floor(Int, I_uout_max - width_fac * width))
    last_idx  = min(length(Ew), floor(Int, I_uout_max + width_fac * width))
    range = first_idx:last_idx

    # Data for polynomial fitting
    pha_x = w[range]
    pha_y = phase_out[range]

    #
    # Third-order polynomial fit
    #
    p3 = fit(pha_x, pha_y, 3)

    # Polynomials.jl stores coefficients in ascending order:
    # p(x) = c0 + c1*x + c2*x² + c3*x³

    GDD = -coeff(p3, 2) * factorial(2)
    TOD = -coeff(p3, 3) * factorial(3)

    println("Calculated to 3rd order:")
    println("    GDD = $GDD ps²")
    println("    TOD = $TOD ps³")

    #
    # Fourth-order polynomial fit
    #
    p4 = fit(pha_x, pha_y, 4)

    GD  = -coeff(p4, 1) * factorial(1)
    GDD = -coeff(p4, 2) * factorial(2)
    TOD = -coeff(p4, 3) * factorial(3)
    FOD = -coeff(p4, 4) * factorial(4)

    println("Calculated to 4th order:")
    println("    GD  = $GD ps")
    println("    GDD = $GDD ps²")
    println("    TOD = $TOD ps³")
    println("    FOD = $FOD ps⁴")

    #
    # Display fitted phase up to third order
    #
    npha_y = coeff(p4,0) .+
             coeff(p4,1) .* pha_x .+
             coeff(p4,2) .* pha_x.^2 .+
             coeff(p4,3) .* pha_x.^3

    # Optional plotting (requires Plots.jl)
    # using Plots
    # plot(pha_x, pha_y,
    #      label="Original",
    #      xlabel="Frequency",
    #      ylabel="Phase")
    # plot!(pha_x, npha_y,
    #       label="4th-order fit (displayed to 3rd)")
    # title!("Fitted to 4th, Displayed to 3rd")

    return phase_out, range, Ew

end