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

using FFTW

"""
    filter_gauss(ui, f3dB, fc, n, fo, df)

Apply an n-th order Gaussian filter in the frequency domain.

# Arguments
- `ui`   : input field amplitude
- `f3dB` : 3 dB bandwidth (THz)
- `fc`   : filter center frequency (THz)
- `n`    : Gaussian filter order
- `fo`   : pulse center frequency (THz)
- `df`   : frequency spacing (THz)

# Returns
- `uo` : filtered field amplitude
"""
function filter_gauss(ui, f3dB, fc, n, fo, df)

    Ui = fft(ui)
    N = length(Ui)

    # Frequency vector (THz)
    f = fftshift(((-N÷2):(N÷2-1)) .* df .+ fo)

    # n-th order Gaussian transfer function
    Tf = @. exp(-log(sqrt(2)) * (2 * (f - fc) / f3dB)^(2n))

    # Apply filter
    uo = ifft(Ui .* Tf)

    return uo
end


"""
    filter_lorentz_tf(ui, fbw, fc, fo, df)

Compute the transfer function of a Lorentzian filter.

# Arguments
- `ui`  : input field amplitude
- `fbw` : Lorentzian filter FWHM (THz)
- `fc`  : filter center frequency (THz)
- `fo`  : pulse center frequency (THz)
- `df`  : frequency spacing (THz)

# Returns
- `tf` : normalized Lorentzian transfer function
"""
function filter_lorentz_tf(ui, fbw, fc, fo, df)

    N = length(ui)

    # Frequency vector (THz)
    f = ((-N÷2):(N÷2-1)) .* df .+ fo

    # Lorentzian transfer function
    tf = @. (fbw / (2π)) / ((f - fc)^2 + (fbw / 2)^2)

    # Normalize to a maximum value of 1
    tf ./= maximum(tf)

    return tf
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
