function calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=false, do_tip_vortex=false, blade_tip=nothing, do_tebvs=false, h=nothing, Psi=nothing)
    if do_tebvs
        if do_tip_vortex
            combined_calc = :with_tip
        else
            combined_calc = :no_tip
        end
    else
        combined_calc = :none
    end

    M_c = 0.8*M
    # How we're going to set this up:
    #
    #   * Source starts at the origin, moving in the -x direction with velocity U.
    #
    c0 = U/M
    # Hmm... how about the location stuff?
    # I want the source element at the origin, so I guess that's r = 0, and θ doesn't matter.
    r = 0.0
    θ = 0.0
    # Hmm... what about the twist?
    # I want the chord line to be at an angle `alphastar` with the negative-x axis.
    # But how does it work "by default"?
    # Let's look at the doc strings.
    # OK, so, if `twist_about_positive_y` is true, then initially the unit vector pointing from leading edge to trailing edge will be pointed in the negative z direction, then rotated about the positive y axis by an amount φ.
    # I want the blade to be aligned with the x axis, with the chord unit vector in the positive x axis direction.
    twist_about_positive_y = true
    # So then I think I want to rotate it by this much:
    ϕ = 3*pi/2 + alphastar
    # ϕ = 3*pi/2
    # OK, so now, how about the velocities?
    # Well, we're starting out in the blade fixed-frame, so the velocity should include everything, including induction.
    # So, from the perspective of the source element, the total velocity is `M_c*c0` in the positive x direction.
    # No velocity in the other directions.
    vn = M_c*c0
    vr = 0.0
    vc = 0.0
    # Ah, now I've decided that's not right.
    # I want the chord line to be aligned with the x axis, so I'll need to rotate the velocity to take into account the angle of attack.
    # So that means the normal velocity would be
    # vn = M_c*c0*cos(alphastar)
    # vr = 0.0
    # vc = M_c*c0*sin(alphastar)
    # We're starting at τ = 0, and the time step doesn't matter yet.
    τ = 0.0
    Δτ = 1.0
    se_tblte = AcousticAnalogies.TBLTESourceElement{AcousticAnalogies.BPMDirectivity,false,AcousticAnalogies.NoMachCorrection}(c0, nu, r, θ, L, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    if do_lblvs
        se_lblvs = AcousticAnalogies.LBLVSSourceElement{AcousticAnalogies.BPMDirectivity,false}(c0, nu, r, θ, L, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    end
    if do_tip_vortex
        se_tip = AcousticAnalogies.TipVortexSourceElement{AcousticAnalogies.BPMDirectivity,false}(c0, r, θ, L, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
    end
    if do_tebvs
        se_tebvs = AcousticAnalogies.TEBVSSourceElement{AcousticAnalogies.BPMDirectivity,false}(c0, nu, r, θ, L, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    end


    if combined_calc == :no_tip
        se_combined = AcousticAnalogies.CombinedNoTipBroadbandSourceElement{AcousticAnalogies.BPMDirectivity,false,AcousticAnalogies.NoMachCorrection}(c0, nu, r, θ, L, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    elseif combined_calc == :with_tip
        se_combined = AcousticAnalogies.CombinedWithTipBroadbandSourceElement{AcousticAnalogies.BPMDirectivity,false,AcousticAnalogies.NoMachCorrection}(c0, nu, r, θ, L, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
    end

    # Let's check that things are what we expect.
    @assert se_tblte.y0dot ≈ [0.0, 0.0, 0.0]
    @assert se_tblte.y1dot ≈ [0.0, 0.0, 0.0]
    @assert se_tblte.y1dot_fluid ≈ [M_c*c0, 0.0, 0.0]
    # @assert se_tblte.y1dot_fluid ≈ [M_c*c0*cos(alphastar), 0.0, M_c*c0*sin(alphastar)]
    @assert se_tblte.span_uvec ≈ [0.0, 1.0, 0.0]
    @assert se_tblte.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
    # @assert se_tblte.chord_uvec ≈ [1.0, 0.0, 0.0]
    @assert isapprox(AcousticAnalogies.angle_of_attack(se_tblte), alphastar; atol=1e-12)

    if do_lblvs
        # Let's check that things are what we expect.
        @assert se_lblvs.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_lblvs.y1dot ≈ [0.0, 0.0, 0.0]
        @assert se_lblvs.y1dot_fluid ≈ [M_c*c0, 0.0, 0.0]
        # @assert se_lblvs.y1dot_fluid ≈ [M_c*c0*cos(alphastar), 0.0, M_c*c0*sin(alphastar)]
        @assert se_lblvs.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_lblvs.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_lblvs.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_lblvs), alphastar; atol=1e-12)
    end

    if do_tip_vortex
        # Let's check that things are what we expect.
        @assert se_tip.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_tip.y1dot ≈ [0.0, 0.0, 0.0]
        @assert se_tip.y1dot_fluid ≈ [M_c*c0, 0.0, 0.0]
        # @assert se_tip.y1dot_fluid ≈ [M_c*c0*cos(alphastar), 0.0, M_c*c0*sin(alphastar)]
        @assert se_tip.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_tip.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_tip.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_tip), alphastar; atol=1e-12)
    end

    if do_tebvs
        # Let's check that things are what we expect.
        @assert se_tebvs.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_tebvs.y1dot ≈ [0.0, 0.0, 0.0]
        @assert se_tebvs.y1dot_fluid ≈ [M_c*c0, 0.0, 0.0]
        # @assert se_tebvs.y1dot_fluid ≈ [M_c*c0*cos(alphastar), 0.0, M_c*c0*sin(alphastar)]
        @assert se_tebvs.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_tebvs.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_tebvs.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_tebvs), alphastar; atol=1e-12)
    end

    if combined_calc != :none
        # Let's check that things are what we expect.
        @assert se_combined.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_combined.y1dot ≈ [0.0, 0.0, 0.0]
        @assert se_combined.y1dot_fluid ≈ [M_c*c0, 0.0, 0.0]
        # @assert se_combined.y1dot_fluid ≈ [M_c*c0*cos(alphastar), 0.0, M_c*c0*sin(alphastar)]
        @assert se_combined.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_combined.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_combined.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_combined), alphastar; atol=1e-12)
    end

    # Now we want to transform the source element into the global frame, which just means we make it move with speed `U` in the negative x direction.
    trans = KinematicCoordinateTransformations.ConstantVelocityTransformation(τ, [0.0, 0.0, 0.0], [-U, 0.0, 0.0])
    # No, that's not right, I want it to move in the direction of the velocity, which isn't aligned with the x axis any more.
    # trans = KinematicCoordinateTransformations.ConstantVelocityTransformation(τ, [0.0, 0.0, 0.0], [-U*cos(alphastar), 0.0, -U*sin(alphastar)])
    se_tblte_global = trans(se_tblte)
    if do_lblvs
        se_lblvs_global = trans(se_lblvs)
    end
    if do_tip_vortex
        se_tip_global = trans(se_tip)
    end
    if do_tebvs
        se_tebvs_global = trans(se_tebvs)
    end
    if combined_calc != :none
        se_combined_global = trans(se_combined)
    end

    # Will that change the angle of attack?
    # Originally the total velocity was `Vtotal = se_tblte.y1dot_fluid - se_tblte.y1dot = [M_c*c0, 0, 0]`, and now it will be `[M_c*c0 - U, 0, 0] - [-U, 0, 0] = [M_c*c0, 0, 0]`.
    # So nope, which is good. :-)
    # Oh, and I'm only changing the magnitude of the velocity, not the direction.
    @assert se_tblte_global.y0dot ≈ [0.0, 0.0, 0.0]
    @assert se_tblte_global.y1dot ≈ [-U, 0.0, 0.0]
    # @assert se_tblte_global.y1dot ≈ [-U*cos(alphastar), 0.0, -U*sin(alphastar)]
    @assert se_tblte_global.y1dot_fluid ≈ [M_c*c0 - U, 0.0, 0.0]
    # @assert se_tblte_global.y1dot_fluid ≈ [(M_c*c0 - U)*cos(alphastar), 0.0, (M_c*c0 - U)*sin(alphastar)]
    @assert se_tblte_global.span_uvec ≈ [0.0, 1.0, 0.0]
    @assert se_tblte_global.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
    # @assert se_tblte_global.chord_uvec ≈ [1.0, 0.0, 0.0]
    @assert isapprox(AcousticAnalogies.angle_of_attack(se_tblte_global), alphastar; atol=1e-12)

    if do_lblvs
        @assert se_lblvs_global.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_lblvs_global.y1dot ≈ [-U, 0.0, 0.0]
        # @assert se_lblvs_global.y1dot ≈ [-U*cos(alphastar), 0.0, -U*sin(alphastar)]
        @assert se_lblvs_global.y1dot_fluid ≈ [M_c*c0 - U, 0.0, 0.0]
        # @assert se_lblvs_global.y1dot_fluid ≈ [(M_c*c0 - U)*cos(alphastar), 0.0, (M_c*c0 - U)*sin(alphastar)]
        @assert se_lblvs_global.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_lblvs_global.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_lblvs_global.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_lblvs_global), alphastar; atol=1e-12)
    end

    if do_tip_vortex
        @assert se_tip_global.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_tip_global.y1dot ≈ [-U, 0.0, 0.0]
        # @assert se_tip_global.y1dot ≈ [-U*cos(alphastar), 0.0, -U*sin(alphastar)]
        @assert se_tip_global.y1dot_fluid ≈ [M_c*c0 - U, 0.0, 0.0]
        # @assert se_tip_global.y1dot_fluid ≈ [(M_c*c0 - U)*cos(alphastar), 0.0, (M_c*c0 - U)*sin(alphastar)]
        @assert se_tip_global.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_tip_global.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_tip_global.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_tip_global), alphastar; atol=1e-12)
    end

    if do_tebvs
        @assert se_tebvs_global.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_tebvs_global.y1dot ≈ [-U, 0.0, 0.0]
        # @assert se_tebvs_global.y1dot ≈ [-U*cos(alphastar), 0.0, -U*sin(alphastar)]
        @assert se_tebvs_global.y1dot_fluid ≈ [M_c*c0 - U, 0.0, 0.0]
        # @assert se_tebvs_global.y1dot_fluid ≈ [(M_c*c0 - U)*cos(alphastar), 0.0, (M_c*c0 - U)*sin(alphastar)]
        @assert se_tebvs_global.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_tebvs_global.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_tebvs_global.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_tebvs_global), alphastar; atol=1e-12)
    end

    if combined_calc != :none
        @assert se_combined_global.y0dot ≈ [0.0, 0.0, 0.0]
        @assert se_combined_global.y1dot ≈ [-U, 0.0, 0.0]
        # @assert se_combined_global.y1dot ≈ [-U*cos(alphastar), 0.0, -U*sin(alphastar)]
        @assert se_combined_global.y1dot_fluid ≈ [M_c*c0 - U, 0.0, 0.0]
        # @assert se_combined_global.y1dot_fluid ≈ [(M_c*c0 - U)*cos(alphastar), 0.0, (M_c*c0 - U)*sin(alphastar)]
        @assert se_combined_global.span_uvec ≈ [0.0, 1.0, 0.0]
        @assert se_combined_global.chord_uvec ≈ [cos(alphastar), 0.0, -sin(alphastar)]
        # @assert se_combined_global.chord_uvec ≈ [1.0, 0.0, 0.0]
        @assert isapprox(AcousticAnalogies.angle_of_attack(se_combined_global), alphastar; atol=1e-12)
    end

    # If the angle of attack is negative, then the pressure and suction sides of the airfoil section switch, and so the coordinate system does too.
    if (alphastar - AcousticAnalogies.alpha_zerolift(bl)) < 0
        Φ_e *= -1
    end

    # What about the observer?
    # That's the tricky part.
    # We know the final position of the observer is this.
    x_obs_final = [r_e*cos(θ_e), r_e*sin(θ_e)*cos(Φ_e), r_e*sin(θ_e)*sin(Φ_e)]
    # This polar coordinate system is actually defined from the perspective of the fluid velocity, not the source element chord direction.
    # So need to take into acount that.
    # x_obs_final = [r_e*cos(θ_e)*cos(alphastar), r_e*sin(θ_e)*cos(Φ_e), r_e*sin(θ_e)*sin(Φ_e)*sin(alphastar)]

    # And we can use the distance from the initial position of the source to the final position of the observer to get the distance the acoustic wave travels, and then the final time.
    t_final = τ + norm(x_obs_final - se_tblte_global.y0dot)/c0

    # And then we can use that to get the initial position of the observer, which is moving with the same velocity as the source.
    x_obs_initial = x_obs_final - se_tblte_global.y1dot*(t_final - τ)

    # So now I should be able to create the observer object thingy.
    obs = AcousticAnalogies.ConstVelocityAcousticObserver(τ, x_obs_initial, se_tblte_global.y1dot)

    # And now I should check if I get the expected final time.
    @assert AcousticAnalogies.adv_time(se_tblte_global, obs) ≈ t_final
    if do_lblvs
        @assert AcousticAnalogies.adv_time(se_lblvs_global, obs) ≈ t_final
    end
    if do_tip_vortex
        @assert AcousticAnalogies.adv_time(se_tip_global, obs) ≈ t_final
    end
    if do_tebvs
        @assert AcousticAnalogies.adv_time(se_tebvs_global, obs) ≈ t_final
    end
    if combined_calc != :none
        @assert AcousticAnalogies.adv_time(se_combined_global, obs) ≈ t_final
    end

    # So now I should be able to do a noise prediction.
    freqs = AcousticMetrics.ExactThirdOctaveCenterBands(0.2, 20e3)
    tblte_out = AcousticAnalogies.noise(se_tblte_global, obs, freqs)
    tblte_s_out = AcousticAnalogies.pbs_suction(tblte_out)
    tblte_p_out = AcousticAnalogies.pbs_pressure(tblte_out)
    tblte_alpha_out = AcousticAnalogies.pbs_alpha(tblte_out)
    SPL_s_jl = 10.0 .* log10.(tblte_s_out./((20e-6)^2))
    SPL_p_jl = 10.0 .* log10.(tblte_p_out./((20e-6)^2))
    SPL_alpha_jl = 10.0 .* log10.(tblte_alpha_out./((20e-6)^2))
    if do_lblvs
        lblvs_out = AcousticAnalogies.noise(se_lblvs_global, obs, freqs)
        SPL_lbl_vs = 10.0 .* log10.(lblvs_out./((20e-6)^2))
    end
    if do_tip_vortex
        tip_out = AcousticAnalogies.noise(se_tip_global, obs, freqs)
        SPL_tip = 10.0 .* log10.(tip_out./((20e-6)^2))
    end
    if do_tebvs
        tebvs_out = AcousticAnalogies.noise(se_tebvs_global, obs, freqs)
        SPL_teb = 10.0 .* log10.(tebvs_out./((20e-6)^2))
    end
    if combined_calc != :none
        combined_out = AcousticAnalogies.noise(se_combined_global, obs, freqs)
    end

    @assert AcousticAnalogies.doppler(tblte_out) ≈ 1
    if do_lblvs
        @assert AcousticAnalogies.doppler(lblvs_out) ≈ 1
    end
    if do_tip_vortex
        @assert AcousticAnalogies.doppler(tip_out) ≈ 1
    end
    if do_tebvs
        @assert AcousticAnalogies.doppler(tebvs_out) ≈ 1
    end
    if combined_calc != :none
        @assert AcousticAnalogies.doppler(combined_out) ≈ 1
    end

    if combined_calc in (:no_tip, :with_tip)
        @assert all(pbs_suction(combined_out) .≈ tblte_s_out)
        @assert all(pbs_pressure(combined_out) .≈ tblte_p_out)
        @assert all(pbs_alpha(combined_out) .≈ tblte_alpha_out)
        @assert all(pbs_teb(combined_out) .≈ tebvs_out)
    end
    if combined_calc in (:with_tip,)
        @assert all(pbs_tip(combined_out) .≈ tip_out)
    end

    res = (freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl)
    if do_lblvs
        res = (res..., SPL_lbl_vs)
    end
    if do_tip_vortex
        res = (res..., SPL_tip)
    end
    if do_tebvs
        res = (res..., SPL_teb)
    end
    return res
end

