function K_1(Re_c)
    # Equation (47) in the BPM paper.
    if Re_c < 2.47e5
        return -4.31*log10(Re_c) + 156.3
    elseif Re_c < 8.0e5
        return -9.0*log10(Re_c) + 181.6
    else
        return 128.5*one(Re_c)
    end
end

function DeltaK_1(alphastar, Re_deltastar_p)
    # Equation (48) in the BPM report.
    T = promote_type(typeof(alphastar), typeof(Re_deltastar_p))
    alphastar_deg = alphastar*180/pi
    if Re_deltastar_p < 5000
        return alphastar_deg*(1.43*log10(Re_deltastar_p) - 5.29)
    else
        return zero(T)
    end
end

function St_1(M)
    # Equation (31) from the BPM report.
    return 0.02*M^(-0.6)
end

function A_min(a)
    # Equation (35) from the BPM report.
    if a < 0.204
        return sqrt(67.552 - 886.788*a^2) - 8.219
    elseif a ≤ 0.244
        return -32.665*a + 3.981
    else
        return -142.795*a^3 + 103.656*a^2 - 57.757*a + 6.006
    end
end

function A_max(a)
    # Equation (36) from the BPM report.
    if a < 0.13
        return sqrt(67.552 - 886.788*a^2) - 8.219
    elseif a ≤ 0.321
        return -15.901*a + 1.098
    else
        return -4.669*a^3 + 3.491*a^2 - 16.699*a + 1.149
    end
end

function A(St_over_St_peak, Re_c)
    # Equation (37) from the BPM report.
    a = abs_cs_safe(log10(St_over_St_peak))
    
    # Equation (38) from the BPM report.
    if Re_c < 9.52e4
        a0 = 0.57*one(Re_c)
    elseif Re_c ≤ 8.57e5
        a0 = (-9.57e-13)*(Re_c - 8.57e5)^2 + 1.13
    else
        a0 = 1.13*one(Re_c)
    end

    # Equation (39) from the BPM report.
    A_min_a0 = A_min(a0)
    A_max_a0 = A_max(a0)
    A_R = (-20 - A_min_a0)/(A_max_a0 - A_min_a0)

    # Equation (40) from the BPM report.
    A_min_a = A_min(a)
    A_max_a = A_max(a)
    return A_min_a + A_R*(A_max_a - A_min_a)
end

function B_min(b)
    # Equation (41) from the BPM report.
    if b < 0.13
        return sqrt(16.888 - 886.788*b^2) - 4.109
    elseif b ≤ 0.145
        return -83.607*b + 8.138
    else
        return -817.810*b^3 + 355.201*b^2 - 135.024*b + 10.619
    end
end

function B_max(b)
    # Equation (42) from the BPM report.
    if b < 0.10
        return sqrt(16.888 - 886.788*b^2) - 4.109
    elseif b ≤ 0.187
        return -31.330*b + 1.854
    else
        return -80.541*b^3 + 44.174*b^2 - 39.381*b + 2.344
    end
end

function B(St_over_St_peak, Re_c)
    # Equation (43) from the BPM report.
    b = abs_cs_safe(log10(St_over_St_peak))

    # Equation (44) from the BPM report.
    if Re_c < 9.52e4
        b0 = 0.30*one(Re_c)
    elseif Re_c ≤ 8.57e5
        b0 = (-4.48e-13)*(Re_c - 8.57e5)^2 + 0.56
    else
        b0 = 0.56*one(Re_c)
    end

    # Equation (45) from the BPM report.
    B_min_b0 = B_min(b0)
    B_max_b0 = B_max(b0)
    B_R = (-20 - B_min_b0)/(B_max_b0 - B_min_b0)

    # Equation (46) from the BPM report.
    B_min_b = B_min(b)
    B_max_b = B_max(b)
    return B_min_b + B_R*(B_max_b - B_min_b)
end

function St_2(St_1, alphastar)
    # Equation (34) from the BPM report.
    T = promote_type(typeof(St_1), typeof(alphastar))
    alphastar_deg = alphastar*180/pi
    if alphastar_deg < 1.333
        return St_1*one(T)
    elseif alphastar_deg ≤ 12.5
        return St_1*10.0^(0.0054*(alphastar_deg - 1.333)^2)
    else
        return St_1*4.72*one(T)
    end
end

function gamma(M)
    # Equation (50) from the BPM report.
    gamma_deg = 27.094*M + 3.31
    return gamma_deg
end

function gamma0(M)
    # Equation (50) from the BPM report.
    gamma0_deg = 23.43*M + 4.651
    return gamma0_deg
end

function beta(M)
    # Equation (50) from the BPM report.
    beta_deg = 72.65*M + 10.74
    return beta_deg
end

function beta0(M)
    # Equation (50) from the BPM report.
    beta0_deg = -34.19*M - 13.82
    return beta0_deg
end

function K_2(Re_c, M, alphastar)
    T = promote_type(typeof(Re_c), typeof(M), typeof(alphastar))
    alphastar_deg = alphastar*180/pi

    k_1 = K_1(Re_c)*one(T)
    # Equation (50) from the BPM report.
    # gamma_deg, gamma0_deg, beta_deg, beta0_deg = gammas_betas(M)
    gamma_deg = gamma(M)
    gamma0_deg = gamma0(M)
    beta_deg = beta(M)
    beta0_deg = beta0(M)

    # Equation (49) from the BPM report.
    if alphastar_deg < gamma0_deg - gamma_deg
        return k_1 - 1000
    # elseif alphastar_deg ≤ gamma0_deg + gamma_deg
    elseif alphastar_deg ≤ gamma0_deg + sqrt(-(gamma_deg/beta_deg)^2*(-12 - beta0_deg)^2 + gamma_deg^2)
        return k_1 + sqrt(beta_deg^2 - (beta_deg/gamma_deg)^2*(alphastar_deg - gamma0_deg)^2) + beta0_deg
    else
        return k_1 - 12
    end
end

function St_1prime(Re_c)
    # Equation (55) from the BPM report.
    T = typeof(Re_c)
    if Re_c < 1.3e5
        return 0.18*one(T)
    elseif Re_c ≤ 4.0e5
        return 0.001756*Re_c^(0.3931)
    else
        return 0.28*one(T)
    end
end

function Dbar_h(theta_e, phi_e, M, M_c)
    # Equation (B1) from the BPM report.
    return (2*sin(0.5*theta_e)^2*sin(phi_e)^2)/((1 + M*cos(theta_e))*(1 + (M - M_c)*cos(theta_e))^2)
end

function Dbar_l(theta_e, phi_e, M)
    # Equation (B2) from the BPM report.
    return (sin(theta_e)^2*sin(phi_e)^2)/((1 + M*cos(theta_e))^4)
end

function TBL_TE_s(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    T = promote_type(typeof(freq), typeof(nu), typeof(L), typeof(chord), typeof(U), typeof(M), typeof(M_c), typeof(r_e), typeof(theta_e), typeof(phi_e), typeof(alphastar))

    Re_c = U*chord/nu
    alphastar0 = alpha_stall(bl, Re_c)

    gamma0_deg = gamma0(M)
    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        # SPL_s = -100*one(T)
        G_s = 10^(0.1*(-100))*one(T)
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)
        deltastar_s = disp_thickness_top(bl, Re_c, alphastar)*chord
        St_s = freq*deltastar_s/U

        St_peak_p = St_1(M)
        St_peak_alpha = St_2(St_peak_p, alphastar)
        St_peak_s = 0.5*(St_peak_p + St_peak_alpha)

        A_s = A(St_s/St_peak_s, Re_c)

        k_1 = K_1(Re_c)

        # SPL_s = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + A_s + k_1 - 3
        # Brooks and Burley AIAA 2001-2210 style.
        H_s = 10^(0.1*(A_s + k_1 - 3))
        G_s = (deltastar_s*M^5*L*D)/(r_e^2)*H_s
    end

    SPL_s = 10*log10(G_s)
    return SPL_s
end

function TBL_TE_p(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    T = promote_type(typeof(freq), typeof(nu), typeof(L), typeof(chord), typeof(U), typeof(M), typeof(M_c), typeof(r_e), typeof(theta_e), typeof(phi_e), typeof(alphastar))

    Re_c = U*chord/nu
    alphastar0 = alpha_stall(bl, Re_c)

    gamma0_deg = gamma0(M)

    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        # SPL_p = -100*one(T)
        G_p = 10^(0.1*(-100))*one(T)
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)
        deltastar_p = disp_thickness_bot(bl, Re_c, alphastar)*chord

        k_1 = K_1(Re_c)
        St_p = freq*deltastar_p/U
        St_peak_p = St_1(M)

        A_p = A(St_p/St_peak_p, Re_c)

        Re_deltastar_p = U*deltastar_p/nu
        Δk_1 = DeltaK_1(alphastar, Re_deltastar_p)

        # SPL_p = 10*log10((deltastar_p*M^5*L*D)/(r_e^2)) + A_p + k_1 - 3 + Δk_1
        # Brooks and Burley AIAA 2001-2210 style.
        H_p = 10^(0.1*(A_p + k_1 - 3 + Δk_1))
        G_p = (deltastar_p*M^5*L*D)/(r_e^2)*H_p
    end

    SPL_p = 10*log10(G_p)
    return SPL_p
end

function TBL_TE_alpha(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    Re_c = U*chord/nu
    alphastar0 = alpha_stall(bl, Re_c)
    gamma0_deg = gamma0(M)
    deltastar_s = disp_thickness_top(bl, Re_c, alphastar)*chord
    St_s = freq*deltastar_s/U
    St_peak_p = St_1(M)
    St_peak_alpha = St_2(St_peak_p, alphastar)
    k2 = K_2(Re_c, M, alphastar)

    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        D = Dbar_l(theta_e, phi_e, M)
        A_prime = A(St_s/St_peak_alpha, 3*Re_c)
        # SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + A_prime + k2
        # Brooks and Burley AIAA 2001-2210 style.
        H_alpha = 10^(0.1*(A_prime + k2))
        G_alpha = (deltastar_s*M^5*L*D)/(r_e^2)*H_alpha
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)
        B_alpha = B(St_s/St_peak_alpha, Re_c)
        # SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + B_alpha + k2
        # Brooks and Burley AIAA 2001-2210 style.
        H_alpha = 10^(0.1*(B_alpha + k2))
        G_alpha = (deltastar_s*M^5*L*D)/(r_e^2)*H_alpha
    end

    SPL_alpha = 10*log10(G_alpha)
    return SPL_alpha
end

function TBL_TE(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    SPL_s = TBL_TE_s(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    SPL_p = TBL_TE_p(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    SPL_alpha = TBL_TE_alpha(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    return SPL_s, SPL_p, SPL_alpha
end

abstract type AbstractDirectivity end
struct BPMDirectivity <: AbstractDirectivity end
struct BrooksBurleyDirectivity <: AbstractDirectivity end

abstract type AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection} <: AbstractCompactSourceElement end

orientation(se::AbstractBroadbandSourceElement) = se.span_uvec

"""
    doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver, t_obs)

Calculate the Doppler shift factor for noise emitted by source element `se` and recieved by observer `obs` at time `t_obs`, i.e. the ratio between an observer frequency `f` and emitted frequency `f_0`.

The correct value for `t_obs` can be found using [`adv_time`](@ref).
"""
function doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver, t_obs)
    # Location of the observer at the observer time.
    x_obs = obs(t_obs)

    # Also need the speed of sound.
    c = se.c0

    # Get a unit vector pointing from the source position at the source time to the observer position at the observer time.
    rv = x_obs .- se.y0dot
    r = norm_cs_safe(rv)
    rhat = rv/r

    # So, now, if I dot the source velocity with `rhat`, that would give me the component of velocity of the source in the direction of the observer, positive if moving toward it, negative if moving away.
    v_src = dot_cs_safe(velocity(se), rhat)

    # And, if I dot the observer velocity `rhat`, that will give me the component of velocity of the observer in the direction of the source, positive if moving *away* from it, negative if moving toward.
    v_obs = dot_cs_safe(velocity(t_obs, obs), rhat)

    # Now we can get the factor.
    factor = (1 - v_obs/c) / (1 - v_src/c)

    return factor
end

"""
    doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver)

Calculate the Doppler shift factor for noise emitted by source element `se` and recieved by observer `obs`, i.e. the ratio between an observer frequency `f` and emitted frequency `f_0`.

The correct value for `t_obs` will be found using [`adv_time`](@ref) internally.
"""
function doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver)
    # Do the advanced time calculation.
    t_obs = adv_time(se, obs)

    return doppler_factor(se, obs, t_obs)
end

function directivity(se::AbstractBroadbandSourceElement{BrooksBurleyDirectivity}, x_obs, top_is_suction)
    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    # But let's do it right, anyway!
    # if se.chord_cross_span_to_get_top_uvec
    #     # But, if the angle of attack is negative, then the "top" of the airfoil (which is normally the suction side) is actually the suction side.
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)
    #     else
    #         z_uvec_tmp = cross(se.span_uvec, se.chord_uvec)
    #     end
    # else
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.span_uvec, se.chord_uvec)
    #     else
    #         z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)
    #     end
    # end
    z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)*ifelse(se.chord_cross_span_to_get_top_uvec, 1, -1)*ifelse(top_is_suction, 1, -1)
    z_uvec = z_uvec_tmp / norm_cs_safe(z_uvec_tmp)
    
    # Component of rv along the chord line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    x_er = dot_cs_safe(rv, se.chord_uvec)

    # Component of rv along the span line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    y_er = dot_cs_safe(rv, se.span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find sin(Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(x_er/r_er)) = sqrt(r_er^2 - x_er^2)/r_er, and so sin(Θ_er)^2 = (r_er^2 - x_er^2)/r_er^2
    sin2Θer = (r_er^2 - x_er^2)/r_er^2

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 + z_er^2)

    # Need to find 2*sin(0.5*Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But there is a half-angle identity that says sin(θ/2)^2 = 0.5*(1 - cos(θ)).
    # So I actually want 2*sin(0.5*Θ_er)^2 = 2*0.5*(1 - cos(Θ_er)) = (1 - cos(Θ_er)).
    # But I can substitute in Θ_er = acos(x_er/r_er) and get 2*sin(0.5*Θ_er)^2 = 1 - x_er/r_er.
    twosin2halfΘer = 1 - x_er/r_er

    # Now just need the denominator: (1 - M_tot*cos(ξR))^4.
    # M_tot is the "total" velocity from... hmm... what perspective?
    # Let's see... it looks like it's suppose to be from the fluid, aka the global frame.
    # The definition is Brooks and Burley AIAA 2001-2210, equation (14):
    #
    #   V_tot = V - V_wt - V_ind
    #
    # where 
    #
    #   * V is the velocity due to the rotation of the blade element
    #   * V_wt is the wind tunnel velocity, which is positive when it goes against the motion of the blade element.
    #   * V_ind is "the induced velocity due to the near and far wake of the rotor," and appears to be positive in roughly the thrust direction.
    #
    # So if I calculate V - V_wt, that's the "actual" velocity of the blade element, i.e., the velocity of the blade element relative to the fluid far away from the blade element, since it doesn't include the induced velocity.
    # That's what I usually think of as the "actual" velocity, since it's what a stationary observer would observe on a calm day.
    # But when we add in the induced velocity, I think what we're finding is the velocity of the blade element relative to the nearfield velocity.
    # Cool.
    # But does that mean I add or subtract `se.y1dot_fluid` from `se.y1dot`?
    # Well, let's think about that.
    # First, let's say I start with stuff in the blade-fixed frame.
    # And let's say I'm imagining that, from the global frame, I'm assuming the
    # blade element is moving in the positive x direction, initially aligned
    # with the y axis
    # So I think all I need to do is just use se.y1dot_fluid + se.y1dot.
    # Now, cos(ξ_r) is defined by equation (18) in Brooks and Burley AIAA 2001-2210, which is the angle between the radiation vector (rv here) and the total velocity (se.y1dot here).
    # But I can simplify that by just finding the unit radiation vector, then dotting that with the velocity vector, and dividing by c0.

    # Unit radiation vector.
    r_uvec = rv./r_er
    # Equation 14 from Brooks and Burley AIAA 2001-2210.
    Vtotal = se.y1dot - se.y1dot_fluid
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Convective amplification factor for the two directivity functions.
    conv_amp = 1/(1 - Mtotcosξr)^4

    # Now I can finally find the directivity function!
    # Equation (19) from Brooks and Burley AIAA 2001-2210.
    # Dl = (sin2Θer*sin2Φer)/(1 - Mtotcosξr)^4
    Dl = (sin2Θer*sin2Φer)*conv_amp

    # Now I can finally find the directivity function!
    # Equation (20) from Brooks and Burley AIAA 2001-2210.
    # Dh = (twosin2halfΘer*sin2Φer)/(1 - Mtotcosξr)^4
    Dh = (twosin2halfΘer*sin2Φer)*conv_amp

    return r_er, Dl, Dh
end

function directivity(se::AbstractBroadbandSourceElement{BPMDirectivity}, x_obs, top_is_suction)
    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # So, the BPM report uses the local flow velocity, not the chord line, to define the x direction.
    # So, I want to get a unit vector in that direction.
    # Should it include induction?
    # It won't matter for comparing to the data in the BPM report, since the flow including and excluding induction would be in the same direction.
    # In the BPM report Appendix B, the x direction is defined as the opposite of the motion of the source element/flat plate, so I guess I won't use induction.
    # But I want the velocity to be normal to the span direction, so let's remove_that.
    # So, want the x direction to be opposite the velocity of the source element.
    x_vec_tmp1 = -se.y1dot
    # Then we want to remove any part of the velocity in the direction of the span.
    x_vec_tmp2 = x_vec_tmp1 - dot_cs_safe(x_vec_tmp1, se.span_uvec)*x_vec_tmp1
    # Now make it a unit vector:
    x_uvec = x_vec_tmp2 / norm_cs_safe(x_vec_tmp2)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    # But it's supposed to be pointing from the pressure to the suction side, which we can figure out, so let's do it the right way.
    # if se.chord_cross_span_to_get_top_uvec
    #     if top_is_suction
    #         z_uvec_tmp = cross(x_uvec, se.span_uvec)
    #     else
    #         z_uvec_tmp = cross(se.span_uvec, x_uvec)
    #     end
    # else
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.span_uvec, x_uvec)
    #     else
    #         z_uvec_tmp = cross(x_uvec, se.span_uvec)
    #     end
    # end
    z_uvec_tmp = cross(x_uvec, se.span_uvec)*ifelse(se.chord_cross_span_to_get_top_uvec, 1, -1)*ifelse(top_is_suction, 1, -1)
    z_uvec = z_uvec_tmp / norm_cs_safe(z_uvec_tmp)
    
    # Component of rv along the chord line (see Figure B3 in the BPM report).
    x_er = dot_cs_safe(rv, x_uvec)

    # Component of rv along the span line (see Figure B3 the BPM report).
    y_er = dot_cs_safe(rv, se.span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find sin(Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(x_er/r_er)) = sqrt(r_er^2 - x_er^2)/r_er, and so sin(Θ_er)^2 = (r_er^2 - x_er^2)/r_er^2
    sin2Θer = (r_er^2 - x_er^2)/r_er^2

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 + z_er^2)

    # Need to find 2*sin(0.5*Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But there is a half-angle identity that says sin(θ/2)^2 = 0.5*(1 - cos(θ)).
    # So I actually want 2*sin(0.5*Θ_er)^2 = 2*0.5*(1 - cos(Θ_er)) = (1 - cos(Θ_er)).
    # But I can substitute in Θ_er = acos(x_er/r_er) and get 2*sin(0.5*Θ_er)^2 = 1 - x_er/r_er.
    twosin2halfΘer = 1 - x_er/r_er

    # Now just need the denominator: (1 - M_tot*cos(ξR))^4.
    # M_tot is the "total" velocity from... hmm... what perspective?
    # Let's see... it looks like it's suppose to be from the fluid, aka the global frame.
    # The definition is Brooks and Burley AIAA 2001-2210, equation (14):
    #
    #   V_tot = V - V_wt - V_ind
    #
    # where 
    #
    #   * V is the velocity due to the rotation of the blade element
    #   * V_wt is the wind tunnel velocity, which is positive when it goes against the motion of the blade element.
    #   * V_ind is "the induced velocity due to the near and far wake of the rotor," and appears to be positive in roughly the thrust direction.
    #
    # So if I calculate V - V_wt, that's the "actual" velocity of the blade element, i.e., the velocity of the blade element relative to the fluid far away from the blade element, since it doesn't include the induced velocity.
    # That's what I usually think of as the "actual" velocity, since it's what a stationary observer would observe on a calm day.
    # But when we add in the induced velocity, I think what we're finding is the velocity of the blade element relative to the nearfield velocity.
    # Cool.
    # But does that mean I add or subtract `se.y1dot_fluid` from `se.y1dot`?
    # Well, let's think about that.
    # First, let's say I start with stuff in the blade-fixed frame.
    # And let's say I'm imagining that, from the global frame, I'm assuming the
    # blade element is moving in the positive x direction, initially aligned
    # with the y axis
    # So I think all I need to do is just use se.y1dot_fluid + se.y1dot.
    # Now, cos(ξ_r) is defined by equation (18) in Brooks and Burley AIAA 2001-2210, which is the angle between the radiation vector (rv here) and the total velocity (se.y1dot here).
    # But I can simplify that by just finding the unit radiation vector, then dotting that with the velocity vector, and dividing by c0.

    # Unit radiation vector.
    r_uvec = rv./r_er
    # For the BPM directivity function, the velocity/Mach number doesn't include induction.
    Vtotal = se.y1dot
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Convective amplification factor for the low-freqency directivity function.
    conv_amp_l = 1/(1 - Mtotcosξr)^4

    # The BPM high-frequency convective amplification factor is a bit different.
    # It has a factor (M - M_c)*cos(Θ_er), which, in the more general coordinate system of Brooks & Burley would be -(M - M_c)*cos(ξ_r).
    # So, the `M` is the speed of the blade element without induction, and `M_c` is the velocity of the blade element including induction.
    # So `M_c = se.y1dot - se.y1dot_fluid` and then `M - M_c = se.y1dot - (se.y1dot - se.y1dot_fluid) = se.y1dot_fluid.
    # And so what we'd want to do is this:
    M_minus_M_ccosξr = dot_cs_safe(se.y1dot_fluid, r_uvec)/se.c0
    conv_amp_h = 1/((1 - Mtotcosξr)*(1 - M_minus_M_ccosξr)^2)

    # Now I can finally find the directivity function!
    # Equation (B2) from the BPM report.
    Dl = (sin2Θer*sin2Φer)*conv_amp_l

    # Now I can finally find the directivity function!
    # Equation (B1) from the BPM report.
    Dh = (twosin2halfΘer*sin2Φer)*conv_amp_h

    return r_er, Dl, Dh
end

function angle_of_attack(se::AbstractBroadbandSourceElement)
    # Find the total velocity of the fluid from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    # Vtotal = -(se.y1dot - se.y1dot_fluid)
    Vtotal = se.y1dot_fluid - se.y1dot 

    # To get the angle of attack, I need to find the components of the velocity in the chordwise direction, and the direction normal to both the chord and span.
    # So, first need to get a vector normal to both the chord and span, pointing from pressure side to suction side.
    normal_uvec_tmp = ifelse(se.chord_cross_span_to_get_top_uvec,
        cross(se.chord_uvec, se.span_uvec),
        cross(se.span_uvec, se.chord_uvec))
    normal_uvec = normal_uvec_tmp ./ norm_cs_safe(normal_uvec_tmp)

    # Now get the component of velocity in the chord_uvec and normal_uvec directions.
    V_chordwise = dot_cs_safe(Vtotal, se.chord_uvec)
    V_normal = dot_cs_safe(Vtotal, normal_uvec)

    # Now we can find the angle of attack.
    alphastar = atan_cs_safe(V_normal, V_chordwise)
    # alphastar = atan(V_normal, V_chordwise)
    
    return alphastar
end

function speed_normal_to_span(se::AbstractBroadbandSourceElement{TDirect,true}) where {TDirect}
    # Find the total velocity of the fluid including induction, from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    Vtotal = se.y1dot_fluid - se.y1dot
    # Find the component of the velocity in the direction of the span.
    Vspan = dot_cs_safe(Vtotal, se.span_uvec)*se.span_uvec
    # Subtract that from the total velocity to get the velocity normal to the span, then get the norm for the speed normal to span.
    return norm_cs_safe(Vtotal - Vspan)
end

function speed_normal_to_span(se::AbstractBroadbandSourceElement{TDirect,false}) where {TDirect}
    # Find the total velocity of the fluid, not including induction, from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    Vtotal = -se.y1dot
    # Find the component of the velocity in the direction of the span.
    Vspan = dot_cs_safe(Vtotal, se.span_uvec)*se.span_uvec
    # Subtract that from the total velocity to get the velocity normal to the span, then get the norm for the speed normal to span.
    return norm_cs_safe(Vtotal - Vspan)
end

function noise(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
    t_obs = adv_time(se, obs)
    return noise(se, obs, t_obs, freqs)
end

abstract type MachCorrection end
struct NoMachCorrection <: MachCorrection end
struct PrandtlGlauertMachCorrection <: MachCorrection end

@concrete struct TBLTESourceElement{TDirect<:AbstractDirectivity,TUInduction,TMachCorrection} <: AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection}
    # Speed of sound, m/s.
    c0
    # Kinematic viscosity, m^2/s
    nu
    # Radial/spanwise length of element, m.
    Δr
    # chord length of element, m.
    chord
    # Source position, m.
    y0dot
    # Source velocity, m/s.
    y1dot
    # Fluid velocity, m/s.
    y1dot_fluid
    # Source time, s.
    τ
    # Time step size, i.e. the amount of time this source element "exists" at with these properties, s.
    Δτ
    # Radial/spanwise unit vector, aka unit vector aligned with the element's span direction.
    span_uvec
    # Chordwise unit vector, aka unit vector aligned with the element's chord line, pointing from leading edge to trailing edge.
    chord_uvec
    # Boundary layer struct, i.e. an AbstractBoundaryLayer.
    bl
    # `Bool` indicating chord_uvec×span_uvec will give a vector pointing from bottom side (usually pressure side) to top side (usually suction side) if `true`, or the opposite if `false`.
    chord_cross_span_to_get_top_uvec
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), and use the PrandtlGlauertMachCorrection.
function TBLTESourceElement(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE) noise using the BPM/Brooks and Burley method, using position and velocity data expressed in a cylindrical coordinate system.

The `r` and `θ` arguments are used to define the radial and circumferential position of the source element in a cylindrical coordinate system.
Likewise, the `vn`, `vr`, and `vc` arguments are used to define the normal, radial, and circumferential velocity of the fluid (in a reference frame moving with the element) in the same cylindrical coordinate system.
The cylindrical coordinate system is defined as follows:

  * The normal/axial direction is in the positive x axis
  * The circumferential/azimuth angle `θ` is defined such that `θ = 0` means the radial direction is aligned with the positive y axis, and a positive `θ` indicates a right-handed rotation around the positive x axis.

The `twist_about_positive_y` is a `Bool` controling how the `ϕ` argument is handled, which in turn controls the orientation of a unit vector defining `chord_uvec` indicating the orientation of the chord line, from leading edge to trailing edge.
If `twist_about_positive_y` is `true`, `chord_uvec` will initially be pointed in the negative-z direction, and then rotated around the positive y axis by an amount `ϕ` before being rotated by the azimuth angle `θ`.
(This would typcially be appropriate for a source element rotating around the positive x axis.)
If `twist_about_positive_y` is `false`, `chord_uvec` will initially be pointed in the positive-z direction, and then rotated around the negative y axis by an amount `ϕ` before being rotated by the azimuth angle `θ`.
(This would typcially be appropriate for a source element rotating around the negative x axis.)

Note that, for a proper noise prediction, the source element needs to be transformed into the "global" frame, aka, the reference frame of the fluid.
This can be done easily with the transformations provided by the `KinematicCoordinateTransformations` package, or manually by modifying the components of the source element struct.

# Arguments
- c0: Ambient speed of sound (m/s)
- nu: Kinematic viscosity (m^2/s)
- r: radial coordinate of the element in the blade-fixed coordinate system (m)
- θ: angular offest of the element in the blade-fixed coordinate system (rad)
- Δr: length of the element (m)
- chord: chord length of blade element (m)
- ϕ: twist of blade element (rad)
- vn: normal velocity of fluid (m/s)
- vr: radial velocity of fluid (m/s)
- vc: circumferential velocity of the fluid (m/s)
- τ: source time (s)
- Δτ: source time duration (s)
- bl: Boundary layer struct, i.e. an AbstractBoundaryLayer.
- twist_about_positive_y: if `true`, apply twist ϕ about positive y axis, negative y axis otherwise
"""
function TBLTESourceElement{TDirect,TUInduction,TMachCorrection}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) where {TDirect,TUInduction,TMachCorrection}
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    y0dot = @SVector [0, r*cθ, r*sθ]
    T = eltype(y0dot)
    y1dot = @SVector zeros(T, 3)
    y1dot_fluid = @SVector [vn, vr*cθ - vc*sθ, vr*sθ + vc*cθ]
    span_uvec = @SVector [0, cθ, sθ]
    if twist_about_positive_y
        chord_uvec = @SVector [-sϕ, cϕ*sθ, -cϕ*cθ]
    else
        chord_uvec = @SVector [-sϕ, -cϕ*sθ, cϕ*cθ]
    end

    chord_cross_span_to_get_top_uvec = twist_about_positive_y
    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), and use the PrandtlGlauertMachCorrection.
function TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::TBLTESourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::TBLTESourceElement{TDirect,TUInduction,TMachCorrection}) where {TDirect,TUInduction,TMachCorrection}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection}(se.c0, se.nu, se.Δr, se.chord, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.chord_cross_span_to_get_top_uvec)
end

doppler(pbs::AcousticMetrics.AbstractProportionalBandSpectrum) = AcousticMetrics.freq_scaler(pbs)

"""
Output of the turbulent boundary layer-trailing edge (TBL-TE) calculation: the acoustic pressure autospectrum centered at time `t` over observer duration `dt` and observer frequencies `cbands` for the suction side `G_s`, pressure side `G_p`, and the separation `G_alpha`.
"""
struct TBLTEOutput{NO,TF,TG<:AbstractVector{TF},TFreqs<:AcousticMetrics.AbstractProportionalBands{NO,:center},TDTime,TTime} <: AcousticMetrics.AbstractProportionalBandSpectrum{NO,TF}
    G_s::TG
    G_p::TG
    G_alpha::TG
    cbands::TFreqs
    dt::TDTime
    t::TTime

    function TBLTEOutput(G_s::TG, G_p::TG, G_alpha::TG, cbands::AcousticMetrics.AbstractProportionalBands{NO,:center}, dt, t) where {NO,TG}
        ncbands = length(cbands)
        length(G_s) == ncbands || throw(ArgumentError("length(G_s) must match length(cbands)"))
        length(G_p) == ncbands || throw(ArgumentError("length(G_p) must match length(cbands)"))
        length(G_alpha) == ncbands || throw(ArgumentError("length(G_alpha) must match length(cbands)"))
        dt > zero(dt) || throw(ArgumentError("dt must be positive"))
        return new{NO,eltype(TG),TG,typeof(cbands),typeof(dt),typeof(t)}(G_s, G_p, G_alpha, cbands, dt, t)
    end
end

@inline function Base.getindex(pbs::TBLTEOutput, i::Int)
    @boundscheck checkbounds(pbs, i)
    return @inbounds pbs.G_s[i] + pbs.G_p[i] + pbs.G_alpha[i]
end

@inline AcousticMetrics.has_observer_time(pbs::TBLTEOutput) = true
@inline AcousticMetrics.observer_time(pbs::TBLTEOutput) = pbs.t
@inline AcousticMetrics.timestep(pbs::TBLTEOutput) = pbs.dt
@inline AcousticMetrics.time_scaler(pbs::TBLTEOutput, period) = timestep(pbs)/period

# function pbs_suction(pbs::TBLTEOutput)
#     t = AcousticMetrics.observer_time(pbs)
#     dt = AcousticMetrics.timestep(pbs)
#     cbands = AcousticMetrics.center_bands(pbs)
#     return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_s, cbands, dt, t)
# end

# function pbs_pressure(pbs::TBLTEOutput)
#     t = AcousticMetrics.observer_time(pbs)
#     dt = AcousticMetrics.timestep(pbs)
#     cbands = AcousticMetrics.center_bands(pbs)
#     return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_p, cbands, dt, t)
# end

# function pbs_alpha(pbs::TBLTEOutput)
#     t = AcousticMetrics.observer_time(pbs)
#     dt = AcousticMetrics.timestep(pbs)
#     cbands = AcousticMetrics.center_bands(pbs)
#     return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_alpha, cbands, dt, t)
# end

function _tble_te_s(freq, deltastar_s_U, Re_c, St_peak_s, k_1, scaler, deep_stall)
    St_s = freq*deltastar_s_U
    A_s = A(St_s/St_peak_s, Re_c)

    # SPL_s = 10*log10((deltastar_s*M^5*L*Dh)/(r_er^2)) + A_s + k_1 - 3
    # Brooks and Burley AIAA 2001-2210 style.
    H_s = 10^(0.1*(A_s + k_1 - 3))
    # G_s = (deltastar_s*M^5*Δr*Dh)/(r_er^2)*H_s
    G_s = scaler*H_s

    return ifelse(deep_stall, 10^(0.1*(-100))*one(typeof(G_s)), G_s)
end

function _tble_te_p(freq, deltastar_p_U, Re_c, St_peak_p, k_1, Δk_1, scaler, deep_stall)

    St_p = freq*deltastar_p_U

    A_p = A(St_p/St_peak_p, Re_c)

    # SPL_p = 10*log10((deltastar_p*M^5*L*Dh)/(r_er^2)) + A_p + k_1 - 3 + Δk_1
    # Brooks and Burley AIAA 2001-2210 style.
    H_p = 10^(0.1*(A_p + k_1 - 3 + Δk_1))
    # G_p = (deltastar_p*M^5*Δr*Dh)/(r_er^2)*H_p
    G_p = scaler*H_p

    return ifelse(deep_stall, 10^(0.1*(-100))*one(typeof(G_p)), G_p)
end

function _tble_te_alpha(freq, Re_c, deltastar_s_U, St_peak_alpha, k_2, scaler_l, scaler_h, deep_stall)
    # Don't know if this is really necessary.
    T = promote_type(typeof(freq), typeof(Re_c), typeof(deltastar_s_U), typeof(St_peak_alpha), typeof(k_2), typeof(scaler_l), typeof(scaler_h))
    St_s = freq*deltastar_s_U

    A_prime_stall = A(St_s/St_peak_alpha, 3*Re_c)
    # SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_er^2)) + A_prime + k_2
    # Brooks and Burley AIAA 2001-2210 style.
    H_alpha_stall = 10^(0.1*(A_prime_stall + k_2))
    # G_alpha_stall = (deltastar_s*M^5*Δr*Dl)/(r_er^2)*H_alpha_stall
    G_alpha_stall = scaler_l*H_alpha_stall*one(T)

    B_alpha = B(St_s/St_peak_alpha, Re_c)
    # SPL_alpha = 10*log10((deltastar_s*M^5*L*Dh)/(r_er^2)) + B_alpha + k_2
    # Brooks and Burley AIAA 2001-2210 style.
    H_alpha = 10^(0.1*(B_alpha + k_2))
    # G_alpha = (deltastar_s*M^5*Δr*Dh)/(r_er^2)*H_alpha
    G_alpha = scaler_h*H_alpha*one(T)

    return ifelse(deep_stall, G_alpha_stall, G_alpha)
end

# Should use traits or something for this.
function mach_correction(se::AbstractBroadbandSourceElement{TDirect,TUInduction,NoMachCorrection}, M) where {TDirect,TUInduction}
    return one(typeof(M))
end

function mach_correction(se::AbstractBroadbandSourceElement{TDirect,TUInduction,PrandtlGlauertMachCorrection}, M) where {TDirect,TUInduction}
    return 1/(1 - M^2)
end

function noise(se::TBLTESourceElement, obs::AbstractAcousticObserver, t_obs, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
    # Position of the observer:
    x_obs = obs(t_obs)

    # Need the angle of attack.
    alphastar = angle_of_attack(se)

    # Need the directivity functions.
    top_is_suction = is_top_suction(se.bl, alphastar)
    r_er, Dl, Dh = directivity(se, x_obs, top_is_suction)

    # Need the fluid velocity normal to the span.
    # Brooks and Burley 2001 are a bit ambiguous on whether it should include induction, or just the freestream and rotation.
    #
    #   * In the nomenclature section: `U` is "flow speed normal to span (`U_mn` with `mn` suppressed).
    #     So that's one point for "no induction."
    #   * In some discussion after equation (8), "The Mach number, `M = U/c0`, represents that component of velocity `U` normal to the span...".
    #     Hard to say one way or the other.
    #   * In equation (12), `U_mn` is the velocity without induction.
    #     So that's another point for "no induction."
    #   * Equation (14) defines `V_tot` as the velocity including the freestream, rotation, and induction.
    #     And then it defines `U` as the part of `V_tot` normal to the span.
    #     So that's a point for "yes induction."
    #   * In the directivity function definitions in equations (19) and (20), `M_tot` is used in the denominator, which seems to make it clear *that* velocity should include induction, since `V_tot` always includes induction.
    #
    # So, at the moment, the TBLTESourceElement type has a parameter TUInduction which, when true, will include induction in the flow speed normal to the span, and not otherwise.
    U = speed_normal_to_span(se)

    # Reynolds number based on chord and the flow speed normal to span.
    Re_c = U*se.chord/se.nu

    # Also need the displacement thicknesses for the pressure and suction sides.
    deltastar_s = disp_thickness_s(se.bl, Re_c, alphastar)*se.chord
    deltastar_p = disp_thickness_p(se.bl, Re_c, alphastar)*se.chord

    # Now that we've decided on the directivity functions and the displacement thickness, and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphastar_positive = abs_cs_safe(alphastar - alpha_zerolift(se.bl))

    # Mach number of the flow speed normal to span.
    M = U/se.c0

    # This stuff is used to decide if the blade element is stalled or not.
    alphastar0 = alpha_stall(se.bl, Re_c)
    gamma0_deg = gamma0(M)
    deep_stall = (alphastar_positive*180/pi) > min(gamma0_deg, alphastar0*180/pi)

    St_peak_p = St_1(M)
    St_peak_alpha = St_2(St_peak_p, alphastar_positive)
    St_peak_s = 0.5*(St_peak_p + St_peak_alpha)

    Re_deltastar_p = U*deltastar_p/se.nu
    k_1 = K_1(Re_c)
    k_2 = K_2(Re_c, M, alphastar_positive)
    Δk_1 = DeltaK_1(alphastar_positive, Re_deltastar_p)

    deltastar_s_U = deltastar_s/U
    deltastar_p_U = deltastar_p/U

    # Brooks and Burley 2001 recommend a Prandtl-Glauert style Mach number correction.
    # But whether or not it's included is dependent on the TMachCorrection type parameter for the source element.
    m_corr = mach_correction(se, M)

    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    G_s_scaler = (deltastar_s*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_s = _tble_te_s.(freqs, deltastar_s_U, Re_c, St_peak_s, k_1, G_s_scaler, deep_stall).*pref2

    G_p_scaler = (deltastar_p*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_p = _tble_te_p.(freqs, deltastar_p_U, Re_c, St_peak_p, k_1, Δk_1, G_p_scaler, deep_stall).*pref2

    G_alpha_scaler_l = (deltastar_s*M^5*se.Δr*Dl)/(r_er^2)*m_corr
    G_alpha_scaler_h = G_s_scaler
    G_alpha = _tble_te_alpha.(freqs, Re_c, deltastar_s_U, St_peak_alpha, k_2, G_alpha_scaler_l, G_alpha_scaler_h, deep_stall).*pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)

    # All done.
    return TBLTEOutput(G_s, G_p, G_alpha, freqs_obs, dt, t_obs)
end
