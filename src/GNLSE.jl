module GNLSE

import Main.Fibra.fiber_types: smf, nasmf, asmf

using FFTW
using Printf

#show(GNLSE)

"""
    IP_CQEM_FD(u0, dt, dz, mod, fo, tol, dplot, quiet)

Solves the Generalized Nonlinear Schrodinger Equation with the complete
Raman response for pulse propagation in an optical fiber using the
Interaction Picture Method combined with the Conserved Quantity Error
method for step-size determination and frequency-domain integration of
the nonlinear operator.

# Arguments
- `u0`    : starting field amplitude (Complex vector)
- `dt`    : time step [ps]
- `dz`    : initial step size
- `mod`   : propagation module parameters (Dict), with keys
            :length      - propagation distance
            :alpha  - power loss coefficient, i.e. P = P0*exp(-alpha*z)
            :gamma  - nonlinearity coefficient
            :beta  - dispersion Taylor coefficients [beta_0 ... beta_m]
- `fo`    : central frequency of the simulation (THz)
- `tol`   : relative photon error
- `dplot` : if 1, plot data will be packed
- `quiet` : if true, suppress progress printout

# Returns
- `u1`       : field at the output
- `nf`       : number of FFTs performed
- `Plotdata` : Dict with saved propagation data (or 0 if dplot != 1)
"""

"""
Need more rigorous documentation include where do we get the values
Also function is not tested if Raman is true

# Arguments
- `t`     : time domain vector
- `mod`   : sinlge mode fiber

# Returns
- `uo` : filtered field amplitude
"""
function IP_CQEM_FD(u0::Vector{ComplexF64}, dt::Float64, dz::Float64, mod::smf, fo::Float64, tol::Float64, dplot::Bool, quiet::Bool)

    println("Starting IP_CQEM_FD")
    nt = length(u0)                                      # number of sample points
    w = fftshift(2*pi .* (-(nt÷2):(nt÷2 - 1)) ./ (dt*nt)) # angular frequencies
    t = collect((-(nt÷2):1:(nt÷2 - 1)) .* dt)                      # time vector (ps)

    # calculate the raman response function in frequency domain
    hrw, fr = Raman_response_w(t::Vector{Float64}, mod::smf)

 # preparation before the IPM
    ufft = fft(u0)
    propagedlength = 0.0
    u1 = copy(u0)
    nf = 1

    if isa(mod, asmf)
        @show(filter)
        gain_w = filter_lorentz_tf(u1, mod.fbw, mod.fc, fo, 1/(dt*nt))
        alpha_0 = mod.alpha
    end

    if dplot == true
        z_all = Float64[]
        ufft_z = Matrix{ComplexF64}(undef, 0, nt)
        u_z = Matrix{ComplexF64}(undef, 0, nt)
    end

    if !quiet
        print("\nSimulation running...      ")
    end

    while propagedlength < mod.length

        if (dz + propagedlength) > mod.length
            dz = mod.length - propagedlength
        end

        # (re)constructing linear operator
        if isa(mod, asmf)
            # modify alpha to include gain
            Pin0 = sum(u1 .* conj.(u1)) / nt
            gain = gain_saturated2(Pin0, mod.gssdB, mod.PsatdBm) .* gain_w
            mod.alpha = alpha_0 .- gain
        end
        
        LOP = Linearoperator_w(mod.alpha, mod.beta, w)


        
    end

    return u1
end

"""
Computes the Raman response (needs to be tested)

# Arguments
- `t`  : time vector
- `mod`  : single mode fiber

# Returns
- `hrw`
- `fr`
"""
function Raman_response_w(t::Vector{Float64}, mod::smf)

    if ~mod.Raman
        return 0.0, 0.0
    else
        # Raman parameters (ps)
        t1 = 12.2e-3
        t2 = 32e-3
        tb = 96e-3
        fc = 0.04
        fb = 0.21
        fa = 1.0 - fc - fb
        fr = 0.245      # fractional contribution of the delayed Raman response
        
        # Time vector beginning at zero
        tres = t .- first(t)
        
        # Raman response components
        ha = ((t1^2 + t2^2) / (t1 * t2^2)) .* exp.(-tres ./ t2) .* sin.(tres ./ t1)
        hb = ((2tb .- tres) ./ tb^2) .* exp.(-tres ./ tb)
        hr = (fa + fc) .* ha .+ fb .* hb
        
        hrw = fft(hr)
        return hrw, fr
    end
end

"""
    filter_lorentz_t(ui::Complex{Float64}, gain_fbw, gain_fc, f0, df)

Compute the transfer function of a Lorentzian filter on the time domain.

# Arguments
- `ui`          : input field complex amplitude
- `gain_fbw`    : gain bandwidth
- `gain_fc`     : gain center freq.
- f0            : center freq. of pulse
- df            : freq. step size 

# Returns
- `tf` : normalized Lorentzian transfer function
"""
function filter_lorentz_tf(ui::Complex{Float64}, gain_fbw::Float64, gain_fc::Float64, f0::Float64, df::Float64)

    N = length(ui)

    # Frequency vector (THz)
    f = ((-N÷2):(N÷2-1)) * df + f0

    # Lorentzian transfer function
    tf = (gain_fbw / (2π)) / ((f - gain_fc)^2 + (gain_fbw/ 2)^2)

    # Normalize to a maximum value of 1
    tf ./= maximum(tf)
    uo = ui .* tf

    return uo
end

""" 
    gain_saturated2(Pin0, gssdB, PsatdBm)
    calculate the gain coefficient of the amplifier given the input power, 
    the small signal gain coefficient and saturation power
    # Arguments
    - Pin: input average power (W)
    - gssdB: small signal gain coefficient(dB)
    - PsdBm : saturation input power(dBm)

    # Returns
    - gain: gain coefficient of the amplifier
"""
function gain_saturated2(Pin0::Float64, gssdB::Float64, PsatdBm::Float64)
    # convert dB to linear
    gss = 10.0^(gssdB/10.0)
    
    # convert dBm to W
    Psat = (10.0^(PsatdBm/10.0))/1000.0
    # calculate the gain coefficient of the amplifier
    gain = gss / (1.0 + Pin0/Psat)
    return gain
end
# fit it belwow
#=
"""
Need to finish this function, it is not tested yet
    gain_saturated3(Pin,GssdB,PsatdBm) 
    calculate the gain of the amplifier given the input power, the small
    signal gain and saturation power 
    # Arguments
    - Pin: input power
    - GssdB: small signal gain(dB)
    - PsatdBm : saturation input power(dBm)
     
    source code by CJH

% G is the saturated gain 
%       G = Gss*exp(-(G-1)Pin/Psat) (eq1)
"""
function gain_saturated3(Pin::Float64, GssdB::Float64, PsatdBm::Float64)

gss = 10.0^(GssdB/10.0)
Psat = (10.0^(PsatdBm/10.0))/1000.0

# numerical calculation of G
G = fzero(@(G)(G-Gss*exp(-(G-1)*Pin/Psat)),Gss/10);

return G
end

"""
Need to finish this function, it is not tested yet
function [ LOP ] = Linearoperator_w( alpha,betaw,w )
%Linearoperator_w Summary of this function goes here
%   Detailed explanation goes here

LOP = -fftshift(alpha/2);
if (length(betaw) == length(w))     % If the user manually specifies beta(w)
    LOP = LOP - 1i*betaw;
    LOP = fftshift(LOP);
else
    for ii = 0:length(betaw)-1;
        LOP = LOP - 1i*betaw(ii+1)*(w).^ii/factorial(ii);
    end
end

end
"""
=#
end