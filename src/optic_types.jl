module  optic_types

include("phys_const.jl")
println("optic_types.jl")

mutable struct 
    lam_c::Float64
    lam_bw::Float64
    lam_fc::Float64
    lamf3dB::Float64
    lam_n::Int
end

function filter_gaussian(lam_c::Float64, lam_bw::Float64, lam_fc::Float64, lam_n::Int)
    return filter_gaussian(lam_c, lam_bw, lam_fc, 0.0, lam_n) 

end

function filter_lorentzian(lam_c::Float64, lam_bw::Float64, lam_fc::Float64, lam_n::Int)
    return filter_lorentzian(lam_c, lam_bw, lam_fc, 0.0, lam_n) 

end

end # module optic_types
