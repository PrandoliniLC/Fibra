



"""
    filter_gauss(ui, f3dB, fc, n, fo, df)

Apply an n-th order Gaussian filter to the input signal.

# Arguments
- `ui`   : Input field amplitude (vector)
- `f3dB` : Gaussian filter 3 dB bandwidth (THz)
- `fc`   : Filter center frequency (THz)
- `n`    : Gaussian filter order
- `fo`   : Pulse center frequency (THz)
- `df`   : Frequency spacing (THz)

# Returns
- `uo` : Filtered output field amplitude
"""
function filter_gauss(ui, f3dB, fc, n, fo, df)

    # FFT of input signal
    Ui = fft(ui)
    N = length(Ui)

    # Frequency vector (THz)
    f = fftshift(((-N÷2):(N÷2-1)) .* df .+ fo)

    # n-th order Gaussian filter transfer function
    Tf = @. exp(-log(sqrt(2)) * (2 * (f - fc) / f3dB)^(2n))

    # Apply filter in frequency domain
    uo = ifft(Ui .* Tf)

    return uo
end

using FFTW

"""
    filter_lorentz_tf(ui, fbw, fc, fo, df)

Compute the transfer function of a Lorentzian filter.

# Arguments
- `ui`  : Input field amplitude (used only to determine the signal length)
- `fbw` : Lorentzian filter bandwidth (FWHM, THz)
- `fc`  : Filter center frequency (THz)
- `fo`  : Pulse center frequency (THz)
- `df`  : Frequency spacing (THz)

# Returns
- `tf` : Normalized Lorentzian transfer function
"""
function filter_lorentz_tf(ui, fbw, fc, fo, df)

    N = length(ui)

    # Frequency vector (THz)
    f = ((-N÷2):(N÷2-1)) .* df .+ fo

    # Lorentzian transfer function
    tf = @. (fbw / (2π)) / ((f - fc)^2 + (fbw / 2)^2)

    # Normalize
    tf ./= maximum(tf)

    return tf
end
"""
    fwhm(x)

Compute the Full Width at Half Maximum (FWHM) of a pulse.

# Arguments
- `x`: Pulse intensity (1D vector)

# Returns
- `width`: FWHM in number of samples
- `I_l`: Left half-maximum index
- `I_r`: Right half-maximum index
"""
function fwhm(x)

    # Peak value and its index
    peak, ind_peak = findmax(x)
    half_peak = peak / 2

    # Left and right portions of the pulse
    x_l = x[1:ind_peak]
    x_r = x[ind_peak+1:end]

    # Distance from the peak to the first half-maximum crossing
    I_l = findfirst(y -> y <= half_peak, reverse(x_l))
    I_r = findfirst(y -> y <= half_peak, x_r)

    # Handle cases where no crossing is found
    isnothing(I_l) && error("Left half-maximum not found.")
    isnothing(I_r) && error("Right half-maximum not found.")

    width = I_l + I_r

    # Convert distances to indices
    I_l = ind_peak - I_l
    I_r = ind_peak + I_r

    return width, I_l, I_r
end