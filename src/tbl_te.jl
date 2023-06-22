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
    a = abs(log10(St_over_St_peak))
    
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
    b = abs(log10(St_over_St_peak))

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
    alphastar0 = stall_alpha(bl, Re_c)

    gamma0_deg = gamma0(M)
    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        # SPL_s = -100*one(T)
        G_s = 10^(0.1*(-100))*one(T)
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)
        deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord
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
    alphastar0 = stall_alpha(bl, Re_c)

    gamma0_deg = gamma0(M)

    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        # SPL_p = -100*one(T)
        G_p = 10^(0.1*(-100))*one(T)
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)
        deltastar_p = disp_thickness_p(bl, Re_c, alphastar)*chord

        k_1 = K_1(Re_c)
        St_p = freq*deltastar_p/U
        St_peak_p = St_1(M)

        A_p = A(St_p/St_peak_p, Re_c)

        Re_deltastar_p = U*deltastar_p/nu
        ΔK_1 = DeltaK_1(alphastar, Re_deltastar_p)

        # SPL_p = 10*log10((deltastar_p*M^5*L*D)/(r_e^2)) + A_p + k_1 - 3 + ΔK_1
        # Brooks and Burley AIAA 2001-2210 style.
        H_p = 10^(0.1*(A_p + k_1 - 3 + ΔK_1))
        G_p = (deltastar_p*M^5*L*D)/(r_e^2)*H_p
    end

    SPL_p = 10*log10(G_p)
    return SPL_p
end

function TBL_TE_alpha(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    Re_c = U*chord/nu
    alphastar0 = stall_alpha(bl, Re_c)
    gamma0_deg = gamma0(M)
    deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord
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

@concrete struct TBLTESourceElement <: AbstractCompactSourceElement
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
    chord_cross_span_to_get_suction_uvec
end

orientation(se::TBLTESourceElement) = se.span_uvec

"""
    TBLTESourceElement(c0, nu, r, θ, Δr, chord, vn, vr, vc, τ, Δτ, bl)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE) noise using the BPM/Brooks and Burley method.

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
function TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y=true)
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    y0dot = @SVector [0, r*cθ, r*sθ]
    T = eltype(y0dot)
    y1dot = @SVector zeros(T, 3)
    y1dot_fluid = @SVector [vn, vr*cθ - vc*sθ, vr*sθ + vc*cθ]
    span_uvec = @SVector [0, cθ, sθ]
    chord_uvec = @SVector [-sϕ, cϕ*sθ, -cϕ*cθ]

    return TBLTESourceElement(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl)
end

function Dbar_l(se::TBLTESourceElement, obs::AcousticObserver, t_obs)
    # Position of the observer.
    x_obs = obs(t_obs)

    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    z_uvec = cross(se.chord_uvec, se.span_uvec)
    
    # Component of rv along the chord line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    x_er = dot_cs_safe(rv, se.chord_uvec)

    # Component of rv along the span line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    y_er = dot_cs_safe(rv, span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find sin(Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(x_er/r_er)) = sqrt(r_er^2 - x_er^2)/r_er, and so sin(Θ_er)^2 = (r_er^2 - x_er^2)/r_er^2
    sin2Θer = (r_er^2 - x_er^2)/r_er^2

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 - z_er^2)

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
    # So I think all I need to do is just use se.y1dot_fluid + se.y1dot.
    # Now, cos(ξ_r) is defined by equation (18) in Brooks and Burley AIAA 2001-2210, which is the angle between the radiation vector (rv here) and the total velocity (se.y1dot here).
    # But I can simplify that by just finding the unit radiation vector, then dotting that with the velocity vector, and dividing by c0.

    # Unit radiation vector.
    r_uvec = rv./r_er
    # Equation 14 from Brooks and Burley AIAA 2001-2210.
    Vtotal = se.y1dot + se.y1dot_fluid
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Now I can finally find the directivity function!
    # Equation (19) from Brooks and Burley AIAA 2001-2210.
    D = (sin2Θer*sin2Φer)/(1 - Mtotcosξr)^4

    return D
end

function Dbar_h(se::TBLTESourceElement, obs::AcousticObserver, t_obs)
    # Position of the observer.
    x_obs = obs(t_obs)

    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    z_uvec = cross(se.chord_uvec, se.span_uvec)
    
    # Component of rv along the chord line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    x_er = dot_cs_safe(rv, se.chord_uvec)

    # Component of rv along the span line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    y_er = dot_cs_safe(rv, span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find 2*sin(0.5*Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But there is a half-angle identity that says sin(θ/2)^2 = 0.5*(1 - cos(θ)).
    # So I actually want 2*sin(0.5*Θ_er)^2 = 2*0.5*(1 - cos(Θ_er)) = (1 - cos(Θ_er)).
    # But I can substitute in Θ_er = acos(x_er/r_er) and get 2*sin(0.5*Θ_er)^2 = 1 - x_er/r_er.
    twosin2halfΘer = 1 - x_er/r_er

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 - z_er^2)

    # Unit radiation vector.
    r_uvec = rv./r_er
    # Equation 14 from Brooks and Burley AIAA 2001-2210.
    Vtotal = se.y1dot + se.y1dot_fluid
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Now I can finally find the directivity function!
    # Equation (20) from Brooks and Burley AIAA 2001-2210.
    D = (twosin2halfΘer*sin2Φer)/(1 - Mtotcosξr)^4

    return D
end

function angle_of_attack(se::TBLTESourceElement)
    # Find the total velocity from the perspective of the blade element, which is just the total velocity of the blade element with the sign switched.
    Vtotal = -(se.y1dot + se.y1dot_fluid)

    # Remove the component of the velocity in the span direction.
    Vtotal_no_span = Vtotal - dot_cs_safe(Vtotal, se.span_uvec)*se.span_uvec

    # Find the angle between Vtotal_no_span and the chord line.
    # The dot product of two vectors is dot(a,b) = |a|*|b|*cos(α).
    # So I can use the dot product to find that angle.
    alphastar = acos(dot_cs_safe(Vtotal_no_span, se.chord_uvec)/norm(Vtotal_no_span))
    
    return alphastar
end

"""
    (trans::KinematicTransformation)(se::TBLTESourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::TBLTESourceElement)
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return TBLTESourceElement(se.c0, se.nu, se.Δr, se.chord, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl)
end

function bpm(se::TBLTESourceElement, obs::AcousticObserver, t_obs)

end

"""
Output of the turbulent boundary layer-trailing edge (TBL-TE) calculation: the acoustic pressure autospectrum centered at time `t` over duration `dt` at frequency `freq` for the suction side `G_s`, pressure side `G_p`, and the separation `G_alpha`.
"""
@concrete struct TBLTEOutput
    t
    dt
    freq
    G_s
    G_p
    G_alpha
end

