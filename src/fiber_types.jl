module fiber_types

include("phys_const.jl")

mutable struct smf # sigle mode fiber
    length::Float64
    alpha::Float64
    gamma::Float64
    beta::Vector{Float64}
    ssp::Bool # supercontinuum generation
    raman::Bool
end 

mutable struct asmf # active single mode fiber
    length::Float64
    alpha::Float64
    gamma::Float64
    beta::Vector{Float64}
    ssp::Bool # supercontinuum generation
    raman::Bool
    # gain-medium fields 
    gssdB::Float64
    PsatdBm::Float64
    lamda_gain::Float64
    lamda_bw::Float64
    fc::Float64
    fbw::Float64
end 

const PM1950 = Dict(
    "CoreNA" => 0.2,                    #     
    "MFD" => 8.0e-6,                    # 1950 nm 
    "Cutoff" => 1720.0e-9,              #     
    "Cladding" => 125.0e-6,             # m
    "CoreD" => 7.0e-6,                  # m
    "keff" => 1.0,                      # 
    "beta" => [0.0, 0.0, -95.0, 0.0],   # 
    "n2" => 30.0,                       # Kerr coefficient in 10^-16 cm^2/W 
    "alpha" => 0.0                      # attenuation in dB/km                         
)

const AMFparam = Dict(
    "small_signal_gain" => 41.0,        # dB              #     
    "saturation_power" => 20.0,         # dBm 
    "gain_bandwidth" => 20.0           # nm [FWHM]      
                          
)

function create_Aeff(fData::Dict, lam::Float64)
    V_parameter = 2*π*(fData["CoreD"]/2)*fData["CoreNA"]/(lam*1e-9)
    MFDc = 2*(fData["CoreD"]/2)*(0.65 + 1.619/V_parameter^(3/2) + 2.879/V_parameter^6)
    return fData["keff"]* π * (MFDc/2)^2
end

function create_smf(fData::Dict, L::Float64, lam::Float64)
    Aeff = create_Aeff(fData, lam)*1e12 # m^2 to um^2
    gamma = 1e4 * 2*π*fData["n2"]/(lam*Aeff) # W^-1 km^-1
    return smf(L, fData["alpha"], gamma, fData["beta"], false, false)
end

function create_asmf(fData::Dict, aData::Dict, L::Float64, lam::Float64)
    Aeff = create_Aeff(fData, lam)*1e12 # m^2 to um^2
    gamma = 1e4 * 2*π*fData["n2"]/(lam*Aeff) # W^-1 km^-1
    return asmf(L, fData["alpha"], gamma, fData["beta"], false, false, aData["small_signal_gain"], aData["saturation_power"], lam, aData["gain_bandwidth"], phys_const.c/lam, phys_const.c*aData["gain_bandwidth"]/lam^2)
end

end # module fiber_types    
