using FFTW
using Printf

"""
    IP_CQEM_FD(u0, dt, dz, mod, fo, tol, dplot, quiet)

Solves the Generalized Nonlinear Schrodinger Equation with the complete
Raman response for pulse propagation in an optical fiber using the
Interaction Picture Method combined with the Conserved Quantity Error
method for step-size determination and frequency-domain integration of
the nonlinear operator.

# Arguments
- `u0`    : starting field amplitude (vector)
- `dt`    : time step [ps]
- `dz`    : initial step size
- `mod`   : propagation module parameters (Dict), with keys
            :L      - propagation distance
            :alpha  - power loss coefficient, i.e. P = P0*exp(-alpha*z)
            :betap  - dispersion Taylor coefficients [beta_0 ... beta_m]
            :betaw  - beta(w) (used internally, built by caller)
            :gamma  - nonlinearity coefficient
- `fo`    : central frequency of the simulation (THz)
- `tol`   : relative photon error
- `dplot` : if 1, plot data will be packed
- `quiet` : if true, suppress progress printout

# Returns
- `u1`       : field at the output
- `nf`       : number of FFTs performed
- `Plotdata` : Dict with saved propagation data (or 0 if dplot != 1)
"""
function IP_CQEM_FD(u0, dt, dz, mod, fo, tol, dplot, quiet)

    nt = length(u0)                                      # number of sample points
    w = fftshift(2*pi .* (-(nt÷2):(nt÷2 - 1)) ./ (dt*nt)) # angular frequencies
    t = (-(nt÷2):1:(nt÷2 - 1)) .* dt                      # time vector (ps)

    # calculate the raman response function in frequency domain
    hrw, fr = Raman_response_w(t, mod)

    # preparation before the IPM
    ufft = fft(u0)
    propagedlength = 0.0
    u1 = copy(u0)
    nf = 1

    if haskey(mod, :gssdB)
        gain_w = filter_lorentz_tf(u1, mod[:fbw], mod[:fc], fo, 1/dt/nt)
        alpha_0 = mod[:alpha]
    end

    if dplot == 1
        z_all = Float64[]
        ufft_z = Vector{Vector{Float64}}()
        u_z = Vector{Vector{ComplexF64}}()
    end

    # Performing the IPM ***********************************************
    if !quiet
        @printf("\nSimulation running...      ")
    end

    while propagedlength < mod[:L]

        if (dz + propagedlength) > mod[:L]
            dz = mod[:L] - propagedlength
        end

        # (re)constructing linear operator
        if haskey(mod, :gssdB)
            # modify alpha to include gain
            Pin0 = sum(u1 .* conj.(u1)) / nt
            gain = gain_saturated2(Pin0, mod[:gssdB], mod[:PsatdBm]) .* gain_w
            mod[:alpha] = alpha_0 .- gain
        end
        LOP = Linearoperator_w(mod[:alpha], mod[:betaw], w)

        # Calculate the real photon number before and after dz propagation
        # neglecting the frequency dependence of S(w) of equation (14) in
        # "Efficient Adaptive Step Size Method for the Simulation of
        # Supercontinuum Generation in Optical Fibers"
        # (Instead of calculating the partial differential of PN in (14),
        # here we directly calculate PN(z))
        PhotonN = sum((abs.(ufft) .^ 2) ./ (w .+ 2*pi*fo))
        PhotonN_z = sum(exp.(-dz .* fftshift(mod[:alpha])) .* (abs.(ufft) .^ 2) ./ (w .+ 2*pi*fo))

        # Applying the Runge-Kutta method in the interaction picture
        # Implements equation (5) in "Optimum Integration Procedures for
        # Supercontinuum Simulation"
        halfstep = exp.(LOP .* dz / 2)
        uip = halfstep .* ufft         # calculate A in interaction picture
        k1 = halfstep .* dz .* NonLinearoperator_w(u1, mod[:gamma], w, fo, fr, hrw, dt, mod)

        uhalf2 = ifft(uip .+ k1 ./ 2)
        k2 = dz .* NonLinearoperator_w(uhalf2, mod[:gamma], w, fo, fr, hrw, dt, mod)

        uhalf3 = ifft(uip .+ k2 ./ 2)
        k3 = dz .* NonLinearoperator_w(uhalf3, mod[:gamma], w, fo, fr, hrw, dt, mod)

        uhalf4 = ifft(halfstep .* (uip .+ k3))
        k4 = dz .* NonLinearoperator_w(uhalf4, mod[:gamma], w, fo, fr, hrw, dt, mod)

        uaux = halfstep .* (uip .+ k1 ./ 6 .+ k2 ./ 3 .+ k3 ./ 3) .+ k4 ./ 6

        propagedlength = propagedlength + dz

        if !quiet
            @printf("\b\b\b\b\b\b%5.2f%%", propagedlength * 100.0 / mod[:L])
        end

        # set dz for the next step
        error = abs(sum((abs.(uaux) .^ 2) ./ (w .+ 2*pi*fo)) - PhotonN_z) / PhotonN_z
        if error > 2 * tol
            # error exceeds double the tolerance, discard this calculation
            # and roll back
            propagedlength = propagedlength - dz
            dz = dz / 2   # reduce the step size by half
        else
            # accept this step and optimize the next step size
            ufft = uaux
            u1 = ifft(ufft)
            if error > tol
                # error exceeds tolerance, reduce step size
                dz = dz / (2^0.2)
            else
                # error is too small, increase step size to accelerate
                # (factor changed from 0.1 to 0.5 to accelerate; speed and
                # precision are now mainly controlled by tolerance)
                if error < 0.5 * tol
                    dz = dz * (2^0.2)
                end
            end
            if dplot == 1
                # save the plot of this step
                push!(z_all, propagedlength)
                push!(ufft_z, abs.(fftshift(ufft)))
                push!(u_z, copy(u1))
            end
        end
        nf = nf + 16
    end

    # pack the output struct
    local Plotdata
    if dplot == 1
        Plotdata = Dict(
            :z    => z_all,
            :ufft => reduce(vcat, transpose.(ufft_z)),
            :u    => abs.(reduce(vcat, transpose.(u_z)))
        )
    else
        Plotdata = 0
    end

    return u1, nf, Plotdata
end