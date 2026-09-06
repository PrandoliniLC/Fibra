module fiber_types

include("phys_const.jl")

abstract type smf end

mutable struct nasmf <: smf # non-active single mode fiber
    length::Float64         # [km] Length
    alpha::Float64          # [1/km] attenuation coefficient
    gamma::Float64          # [1/(W km)] fiber nonlinear coefficient
    beta::Vector{Float64}      # beta coefficients (ps^n/km)
    ssp::Bool # supercontinuum generation
    Raman::Bool
end 

mutable struct asmf <: smf # active single mode fiber
    length::Float64         # [km] Length 
    alpha::Float64          # [1/km] attenuation coefficient
    gamma::Float64          # [1/(W km)] fiber nonlinear coefficient
    beta::Vector{Float64}   # beta coefficients (ps^n/km)
    ssp::Bool               # supercontinuum generation
    Raman::Bool             # Raman active
    # gain-medium fields 
    gssdB::Float64          # [dB] small signal gain coefficient
    PsatdBm::Float64        # [dBm] saturation input power 
    lam_gain::Float64       # [nm] center of gain 
    lam_bw::Float64         # [nm] gain bandwidth
    fc::Float64             # [THz] center of gain 
    fbw::Float64            # [THz] bandwidth of gain 
end 

# https://www.coherent.com/components-accessories/specialty-optical-fibers/single-mode/SM1950
# keff(wavelength) is the mapping value see ECOC '95, paper Tu.L.2.4. Petermann II method
# Also see "Mode Field Diameter and Effective Area" White Paper from Corning  
const PM1950 = Dict(
    "CoreNA" => 0.2,                    #  Numerical Aperture
    "MFD" => 8.0e-6,                    # [m] Gaussian at 1950 nm 
    "Cutoff" => 1720.0e-9,              # [m]    
    "Cladding" => 125.0e-6,             # [m] 
    "CoreD" => 7.0e-6,                  # [m]
    "keff" => 1.0,                      # see above  
    "beta" => [0.0, 0.0, -95.0, 0.0],   # [ps^n/km] [beta0, beta1, beta2, beta3] OL47(4) 822 (2022)   
    "n2" => 30.0,                       # Nonlinear refractive index, [1e-16 cm^2/W] 
    "alpha" => 0.0                      # [1/km]                         
)

# Needs to be checked with the real data from the manufacturer
const AMFparam = Dict(
    "small_signal_gain" => 41.0,        # dB                   
    "saturation_power" => 20.0,         # dBm 
    "gain_bandwidth" => 20.0            # [nm] [FWHM]      
)

"""
Aeff: see "Mode Field Diameter and Effective Area" White Paper from Corning
V-parameter: see https://www.rp-photonics.com/v_number.html
MFDc: Marcuse's equation see https://www.rp-photonics.com/mode_radius.html
"""
function create_Aeff(fData::Dict, lam::Float64)
    V_parameter = 2*π*(fData["CoreD"]/2)*fData["CoreNA"]/(lam*1e-9)
    MFDc = 2*(fData["CoreD"]/2)*(0.65 + 1.619/V_parameter^(3/2) + 2.879/V_parameter^6)
    return fData["keff"]* π * (MFDc/2)^2
end

"""
Create a single-mode fiber (SMF) instance with the specified parameters.
Aeff [µm^2], gamma [1/(W km)], beta [ps^n/km]
L [km]
lam [nm] 
"""
function create_nasmf(fData::Dict, L::Float64, lam::Float64)
    Aeff = create_Aeff(fData, lam)*1e12 # m^2 to um^2
    gamma = 1e4 * 2*π*fData["n2"]/(lam*Aeff) # W^-1 km^-1
    return nasmf(L, fData["alpha"], gamma, fData["beta"], false, true)
end

"""
Create an active single-mode fiber (SMF) instance with the specified parameters.
Aeff [µm^2], gamma [1/(W km)], beta [ps^n/km]
L [km]
lam [nm] 
"""
function create_asmf(fData::Dict, aData::Dict, L::Float64, lam::Float64, lam_gain::Float64)
    Aeff = create_Aeff(fData, lam)*1e12 # m^2 to um^2
    gamma = 1e4 * 2*π*fData["n2"]/(lam*Aeff) # W^-1 km^-1
    return asmf(L, fData["alpha"], gamma, fData["beta"], false, false, aData["small_signal_gain"], aData["saturation_power"], lam_gain, aData["gain_bandwidth"], phys_const.c/lam_gain, phys_const.c*aData["gain_bandwidth"]/lam_gain^2)
end

end # module fiber_types    

"""
    coupler(u1i, u2i, rho)

Directional coupler.

# Arguments
- `u1i`: Complex Input field 1
- `u2i`: Complex Input field 2
- `rho`: Coupling ratio (clamped to [0, 1])

# Returns
- `u1o`: Complex Output field 1
- `u2o`: Complex Output field 2
"""
function coupler(u1i::Complex{Float64}, u2i::Complex{Float64}, rho::Float64)
    # Clamp rho to the interval [0, 1]
    rho = clamp(rho, 0.0, 1.0)

    u1o = sqrt(rho) * u1i + im * sqrt(1 - rho) * u2i
    u2o = im * sqrt(1 - rho) * u1i + sqrt(rho) * u2i

    return u1o, u2o
end