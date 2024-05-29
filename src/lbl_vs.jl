function St_1_prime(Re_c)
    # Equation 55 in the BPM report.
    T = typeof(Re_c)
    if Re_c ≤ 1.3e5
        return 0.18*one(T)
    elseif Re_c ≤ 4.0e5
        return 0.001756*Re_c^(0.3931)
    else
        return 0.28*one(T)
    end
end

function St_peak_prime(St_1_p, alphastar)
    # Equation 56 in the BPM report.
    alphastar_deg = alphastar*180/pi
    return St_1_p*10.0^(-0.04*alphastar_deg)
end

function G1(St_prime_over_St_peak_prime)
    # Equation 57 in the BPM report.
    e = St_prime_over_St_peak_prime
    if e ≤ 0.5974
        return 39.8*log10(e) - 11.12
    elseif e ≤ 0.8545
        98.409*log10(e) + 2.0
    elseif e ≤ 1.17
        return -5.076 + sqrt(2.484 - 506.25*log10(e)^2)
    elseif e ≤ 1.674
        return -98.409*log10(e) + 2.0
    else
        return -39.8*log10(e) - 11.12
    end
end

function Re_c0(alphastar)
    # Equation 59 in the BPM report.
    alphastar_deg = alphastar*180/pi
    if alphastar_deg ≤ 3.0
        return 10.0^(0.215*alphastar_deg + 4.978)
    else
        return 10.0^(0.120*alphastar_deg + 5.263)
    end
end

function G2(Re_c_over_Re_c0)
    # Equation 58 in the BPM report.
    d = Re_c_over_Re_c0
    if d ≤ 0.3237
        return 77.852*log10(d) + 15.328
    elseif d ≤ 0.5689
        return 65.188*log10(d) + 9.125
    elseif d ≤ 1.7579
        return -114.052*log10(d)^2
    elseif d ≤ 3.0889
        return -65.188*log10(d) + 9.125
    else
        return -77.852*log10(d) + 15.328
    end
end

function G3(alphastar)
    alphastar_deg = alphastar*180/pi
    return 171.04 - 3.03*alphastar_deg
end

function LBL_VS(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    Re_c = U*chord/nu
    delta_p = bl_thickness_p(bl, Re_c, alphastar)*chord
    # Equation (54) from the BPM report.
    St_prime = freq*delta_p/U
    St_p_p = St_peak_prime(St_1_prime(Re_c), alphastar)

    St_prime_over_St_peak_prime = St_prime/St_p_p
    Re_c_over_Re_c0 = Re_c / Re_c0(alphastar)
    # SPL = 10*log10(delta_p*M^5*L*Dbar_h(theta_e, phi_e, M, M_c)/r_e^2) + G1(St_prime_over_St_peak_prime) + G2(Re_c_over_Re_c0) + G3(alphastar)
    # Brooks and Burley AIAA 2001-2210 style.
    H_l = 10^(0.1*(G1(St_prime_over_St_peak_prime) + G2(Re_c_over_Re_c0) + G3(alphastar)))
    G_lbl_vs = (delta_p*M^5*L*Dbar_h(theta_e, phi_e, M, M_c))/(r_e^2)*H_l
    SPL = 10*log10(G_lbl_vs)
    return SPL
end

@concrete struct LBLVSSourceElement{TDirect<:AbstractDirectivity,TUInduction} <: AbstractBroadbandSourceElement{TDirect,TUInduction,NoMachCorrection}
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

# Default to using the `BrooksBurleyDirectivity` directivity function, and include induction in the flow speed normal to span (TUInduction == true).
function LBLVSSourceElement(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
    return LBLVSSourceElement{BrooksBurleyDirectivity,true}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    LBLVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting laminar boundary layer-vortex shedding (LBLVS) noise using the BPM/Brooks and Burley method, using position and velocity data expressed in a cylindrical coordinate system.

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
function LBLVSSourceElement{TDirect,TUInduction}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) where {TDirect,TUInduction}
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
    return LBLVSSourceElement{TDirect,TUInduction}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, and include induction in the flow speed normal to span (TUInduction == true).
function LBLVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    return LBLVSSourceElement{BrooksBurleyDirectivity,true}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::LBLVSSourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::LBLVSSourceElement{TDirect,TUInduction}) where {TDirect,TUInduction}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return LBLVSSourceElement{TDirect,TUInduction}(se.c0, se.nu, se.Δr, se.chord, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.chord_cross_span_to_get_top_uvec)
end

function _lbl_vs(freq, delta_p_U, St_p_p, g2, g3, scaler)
    # St_prime = freq*deltastar_p/U
    # St_prime_over_St_peak_prime = St_prime/St_p_p
    # H_l = 10^(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3))
    St_prime = freq*delta_p_U
    St_prime_over_St_peak_prime = St_prime/St_p_p
    # return 10^(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3))
    H_l = 10^(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3))
    # G_lbl_vs = (deltastar_p*M^5*Δr*Dh)/(r_er^2)*H_l
    G_lbl_vs = scaler*H_l

    # LBLVS = 10.0*log10((dpr*M^5*L*Dh)/rc^2)+G1+G2+G3
    # 10*log10((deltastar_p*M^5*Δr*Dh)/(r_er^2)*(10^(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3)))
    # 10*log10((deltastar_p*M^5*Δr*Dh)/(r_er^2)) + 10*log10((10^(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3))))
    # 10*log10((deltastar_p*M^5*Δr*Dh)/(r_er^2)) + 10*(0.1*(G1(St_prime_over_St_peak_prime) + g2 + g3))
    # 10*log10((deltastar_p*M^5*Δr*Dh)/(r_er^2)) + G1(St_prime_over_St_peak_prime) + g2 + g3
    return G_lbl_vs
end

function noise(se::LBLVSSourceElement, obs::AbstractAcousticObserver, t_obs, freqs::AcousticMetrics.ExactThirdOctaveCenterBands)
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

    # Need the boundary layer thickness for the pressure side for LBL-VS noise.
    # deltastar_p = disp_thickness_p(se.bl, Re_c, alphastar)*se.chord
    delta_p = bl_thickness_p(se.bl, Re_c, alphastar)*se.chord

    # Now that we've decided on the directivity functions and the displacement thickness, and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphastar_positive = abs_cs_safe(alphastar - alpha_zerolift(se.bl))

    # Mach number of the flow speed normal to span.
    M = U/se.c0

    delta_p_U = delta_p/U
    St_p_p = St_peak_prime(St_1_prime(Re_c), alphastar_positive)
    Re_c_over_Re_c0 = Re_c / Re_c0(alphastar_positive)
    g2 = G2(Re_c_over_Re_c0)
    g3 = G3(alphastar_positive)
    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    # G = (delta_p*M^5*se.Δr*Dh)/(r_er^2) .* _Hl.(freq, delta_p_U, St_p_p, g2, g3) .* pref2
    G_lbl_vs_scaler = (delta_p*M^5*se.Δr*Dh)/(r_er^2)
    G_lbl_vs = _lbl_vs.(freqs, delta_p_U, St_p_p, g2, g3, G_lbl_vs_scaler) .* pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)

    # All done.
    return AcousticMetrics.ProportionalBandSpectrumWithTime(G_lbl_vs, freqs_obs, dt, t_obs)
end
