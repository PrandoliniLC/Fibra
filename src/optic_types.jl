module  optic_types
import FFTW

mutable struct filter
    f0::Float64             # THz Center Freq. 
    BW::Float64             # [nm] Bandwidth FWHM
    df::Float64             # Thz Frequency spacing
    n::Int                  # Super-Gaussian order n=1,2,3,4
end

"""
    filter_gauss(ui::Complex{Float64}, fl::filter)

Apply an n-th order Gaussian filter in the frequency domain.

# Arguments
- `ui`   : input field amplitude
- `fl`   : filter parameters (type `filter`)

# Returns
- `uo` : filtered field amplitude
"""
function filter_gauss(ui::Complex{Float64}, fl::filter)

    Ui = fft(ui)
    N = length(Ui)

    # Frequency vector (THz)
    f = fftshift(((-N÷2):(N÷2-1)) * fl.df + fl.fo)

    # n-th order Gaussian transfer function
    Tf = exp(-log(sqrt(2)) * (2 * (f - fl.fc) / fl.BW)^(2fl.n))

    # Apply filter
    uo = ifft(Ui .* Tf)

    return uo
end

#=
"""
    filter_lorentz_t(ui::Complex{Float64}, fl::filter)

Compute the transfer function of a Lorentzian filter on the time domain.

# Arguments
- `ui`  : input field amplitude
- `fl`  : filter parameters (type `filter`)

# Returns
- `tf` : normalized Lorentzian transfer function
"""
function filter_lorentz_tf(ui::Complex{Float64}, fl::filter)

    N = length(ui)

    # Frequency vector (THz)
    f = ((-N÷2):(N÷2-1)) * fl.df + fl.fo

    # Lorentzian transfer function
    tf = (fl.BW / (2π)) / ((f - fl.fc)^2 + (fl.BW/ 2)^2)

    # Normalize to a maximum value of 1
    tf ./= maximum(tf)
    uo = ui .* tf

    return uo
end
=#


end # module optic_types