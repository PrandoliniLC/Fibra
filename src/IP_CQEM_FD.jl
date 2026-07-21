function IP_CQEM_FD(
    u0,
    dt,
    dz,
    mod,
    fo,
    tol,
    dplot,
    quiet
)

    # Number of sample points
    nt = length(u0)

    # Angular frequency vector
    w = fftshift(2π .* (-nt÷2:nt÷2-1) ./ (dt * nt))

    # Time vector
    t = (-nt÷2:nt÷2-1) .* dt

    # Raman response
    hrw, fr = Raman_response_w(t, mod)

    # Initialization
    ufft = fft(u0)
    propagated_length = 0.0
    u1 = copy(u0)
    nf = 1

    if hasproperty(mod, :gssdB)
        gain_w = filter_lorentz_tf(
            u1,
            mod.fbw,
            mod.fc,
            fo,
            1 / (dt * nt)
        )
        alpha0 = mod.alpha
    end

    if dplot == 1
        z_all = Float64[]
        ufft_z = Matrix{Float64}(undef, 0, nt)
        u_z = Matrix{ComplexF64}(undef, 0, nt)
    end

    if !quiet
        print("\nSimulation running...      ")
    end

    while propagated_length < mod.L

        if propagated_length + dz > mod.L
            dz = mod.L - propagated_length
        end

        ############################################################
        # Construct linear operator
        ############################################################

        if hasproperty(mod, :gssdB)

            Pin0 = sum(abs2, u1) / nt

            gain = gain_saturated2(
                Pin0,
                mod.gssdB,
                mod.PsatdBm
            ) .* gain_w

            mod.alpha = alpha0 .- gain
        end

        LOP = Linearoperator_w(mod.alpha, mod.betaw, w)

        ############################################################
        # Photon number
        ############################################################

        PhotonN =
            sum(abs2.(ufft) ./ (w .+ 2π * fo))

        PhotonN_z =
            sum(exp.(-dz .* fftshift(mod.alpha))
                .* abs2.(ufft)
                ./ (w .+ 2π * fo))

        ############################################################
        # RK4 in interaction picture
        ############################################################

        halfstep = exp.(LOP .* dz / 2)

        uip = halfstep .* ufft

        k1 = halfstep .* dz .*
             NonLinearoperator_w(
                 u1,
                 mod.gamma,
                 w,
                 fo,
                 fr,
                 hrw,
                 dt,
                 mod
             )

        uhalf2 = ifft(uip .+ k1 ./ 2)

        k2 = dz .* NonLinearoperator_w(
            uhalf2,
            mod.gamma,
            w,
            fo,
            fr,
            hrw,
            dt,
            mod
        )

        uhalf3 = ifft(uip .+ k2 ./ 2)

        k3 = dz .* NonLinearoperator_w(
            uhalf3,
            mod.gamma,
            w,
            fo,
            fr,
            hrw,
            dt,
            mod
        )

        uhalf4 = ifft(halfstep .* (uip .+ k3))

        k4 = dz .* NonLinearoperator_w(
            uhalf4,
            mod.gamma,
            w,
            fo,
            fr,
            hrw,
            dt,
            mod
        )

        uaux =
            halfstep .* (uip .+ k1 ./ 6 .+ k2 ./ 3 .+ k3 ./ 3) .+
            k4 ./ 6

        propagated_length += dz

        if !quiet
            print("\r$(round(propagated_length / mod.L * 100; digits=2)) %")
        end

        ############################################################
        # Adaptive step size
        ############################################################

        err =
            abs(
                sum(abs2.(uaux) ./ (w .+ 2π * fo)) - PhotonN_z
            ) / PhotonN_z

        if err > 2tol

            propagated_length -= dz
            dz /= 2

        else

            ufft = uaux
            u1 = ifft(ufft)

            if err > tol
                dz /= 2^0.2
            elseif err < 0.5tol
                dz *= 2^0.2
            end

            if dplot == 1
                push!(z_all, propagated_length)

                ufft_z = vcat(
                    ufft_z,
                    reshape(abs.(fftshift(ufft)), 1, :)
                )

                u_z = vcat(
                    u_z,
                    reshape(u1, 1, :)
                )
            end
        end

        nf += 16
    end

    ############################################################
    # Output plotting data
    ############################################################

    Plotdata =
        dplot == 1 ?
        (
            z = z_all,
            ufft = ufft_z,
            u = abs.(u_z)
        ) :
        nothing

    return u1, nf, Plotdata

end