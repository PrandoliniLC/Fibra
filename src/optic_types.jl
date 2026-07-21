module  optic_types
using FFTW

include("phys_const.jl")

mutable struct 
    lam_c::Float64
    lam_bw::Float64
    lam_fc::Float64
    lamf3dB::Float64
    lam_n::Int
end

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



end # module optic_types
