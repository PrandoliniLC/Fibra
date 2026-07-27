
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