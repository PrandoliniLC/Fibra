module nonlinear

"""
Need more rigorous documentation include where do we get the values
Also function is not tested if Raman is true

# Arguments
- `t`     : time domain vector
- `mod`   : sinlge mode fiber

# Returns
- `uo` : filtered field amplitude
"""
function Raman_response_w(t::Vector(Float64), mod::smf)

    # Raman response disabled
    if ~mod.Raman
        return 0.0, 0.0
    end

    # Raman parameters (ps)
    t1 = 12.2e-3
    t2 = 32e-3
    tb = 96e-3

    fc = 0.04
    fb = 0.21
    fa = 1.0 - fc - fb

    fr = 0.245

    # Time vector beginning at zero
    tres = t .- first(t)

    # Raman response components
    ha = ((t1^2 + t2^2) / (t1 * t2^2)) .*
         exp.(-tres ./ t2) .*
         sin.(tres ./ t1)

    hb = ((2tb .- tres) ./ tb^2) .*
         exp.(-tres ./ tb)

    hr = (fa + fc) .* ha .+ fb .* hb

    hrw = fft(hr)

    return hrw, fr

end


end