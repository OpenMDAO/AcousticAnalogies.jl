function St_3prime_peak(h_over_deltastar_avg, Psi)
    Psi_deg = Psi*180/pi
    # Equation 72 from the BPM report.
    if h_over_deltastar_avg < 0.2
        return 0.1*(h_over_deltastar_avg) + 0.095 - 0.00243*Psi_deg
    else
        return (0.212 - 0.0045*Psi_deg)/(1 + 0.235/h_over_deltastar_avg - 0.0132/(h_over_deltastar_avg^2))
    end
end

function G4(h_over_deltastar_avg, Psi)
    Psi_deg = Psi*180/pi
    # Equation 74 from the BPM report.
    if h_over_deltastar_avg ≤ 5
        return 17.5*log10(h_over_deltastar_avg) + 157.5 - 1.114*Psi_deg
    else
        return 169.7 - 1.114*Psi_deg
    end
end

function G5_Psi14(h_over_deltastar_avg, St_3prime_over_St_3prime_peak)
    # Equation 77 from the BPM report.
    η = log10(St_3prime_over_St_3prime_peak)

    # Equation 78 from the BPM report.
    if h_over_deltastar_avg < 0.25
        μ = 0.1211*one(h_over_deltastar_avg)
    elseif h_over_deltastar_avg < 0.62
        μ = -0.2175*h_over_deltastar_avg + 0.1755
    elseif h_over_deltastar_avg < 1.15
        μ = -0.0308*h_over_deltastar_avg + 0.0596
    else
        μ = 0.0242*one(h_over_deltastar_avg)
    end

    # Equation 79 from the BPM report.
    if h_over_deltastar_avg < 0.02
        m = zero(h_over_deltastar_avg)
    elseif h_over_deltastar_avg ≤ 0.5
        m = 68.724*h_over_deltastar_avg - 1.35
    elseif h_over_deltastar_avg ≤ 0.62
        m = 308.475*h_over_deltastar_avg - 121.23
    elseif h_over_deltastar_avg < 1.15
        m = 224.811*h_over_deltastar_avg - 69.354
    elseif h_over_deltastar_avg < 1.2
        m = 1583.28*h_over_deltastar_avg - 1631.592
    else
        m = 268.344*one(h_over_deltastar_avg)
    end
    # This is in the code listing in the BPM report appendix.
    if m < 0
        m = zero(h_over_deltastar_avg)
    end

    # Equation 80 from the BPM report.
    η_0 = -sqrt((m^2*μ^4)/(6.25 + m^2*μ^2))

    # Equation 81 from the BPM report.
    k = 2.5*sqrt(1 - (η_0/μ)^2) - 2.5 - m*η_0

    # Equation 76 from the BPM report.
    if η < η_0
        return m*η + k
    elseif η < 0
        return 2.5*sqrt(1 - (η/μ)^2) - 2.5
    elseif η < 0.03616
        return sqrt(1.5625 - 1194.99*η^2) - 1.25*one(h_over_deltastar_avg)
    else
        return -155.543*η + 4.375*one(h_over_deltastar_avg)
    end
end

function G5_Psi0(h_over_deltastar_avg, St_3prime_over_St_3prime_peak)
    # Equation 82 from the BPM report.
    h_over_deltastar_avg_prime = 6.724*h_over_deltastar_avg^2 - 4.019*h_over_deltastar_avg + 1.107
    return G5_Psi14(h_over_deltastar_avg_prime, St_3prime_over_St_3prime_peak)
end

function G5(h_over_deltastar_avg, Psi, St_3prime_over_St_3prime_peak)
    Psi_deg = 180/pi*Psi
    # Equation 75 from the BPM report.
    G5_0 = G5_Psi0(h_over_deltastar_avg, St_3prime_over_St_3prime_peak)
    G5_14 = G5_Psi14(h_over_deltastar_avg, St_3prime_over_St_3prime_peak)
    g5 = G5_0 + 0.0714*Psi_deg*(G5_14 - G5_0)
    # This check is in the code listing in the BPM report appendix:
    if g5 > 0
        return zero(g5)
    else
        return g5
    end
end

function BLUNT(freq, nu, L, chord, h, Psi, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    Re_c = U*chord/nu
    # deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord
    # deltastar_p = disp_thickness_p(bl, Re_c, alphastar)*chord
    deltastar_top = disp_thickness_top(bl, Re_c, alphastar)*chord
    deltastar_bot = disp_thickness_bot(bl, Re_c, alphastar)*chord
    top_is_suction = alphastar > alpha_zerolift(bl)
    deltastar_s, deltastar_p = ifelse(top_is_suction,
        (deltastar_top, deltastar_bot),
        (deltastar_bot, deltastar_top))

    # Equation 73 from the BPM report.
    deltastar_avg = 0.5*(deltastar_p + deltastar_s)

    h_over_deltastar_avg = h/deltastar_avg
    D = Dbar_h(theta_e, phi_e, M, M_c)

    # Equation 71 from the BPM report.
    St_3p = freq*h/U

    St_3prime_over_St_3prime_peak = St_3p/St_3prime_peak(h_over_deltastar_avg, Psi)

    g5temp = G5(h_over_deltastar_avg, Psi, St_3prime_over_St_3prime_peak)

    # This next check is in the code listing in the BPM report appendix.
    # Need to find G5 for h_over_deltastar_avg = 0.25 for the F4TEMP variable.
    f4temp = G5_Psi14(0.25, St_3prime_over_St_3prime_peak)
    # if g5 > f4temp
    #     g5 = f4temp
    # end
    g5 = ifelse(g5temp > f4temp, f4temp, g5temp)

    # Equation 70 from the BPM report.
    # SPL_blunt = 10*log10((h*(M^5.5)*L*D)/(r_e^2)) + G4(h_over_deltastar_avg, Psi) + g5
    # Brooks and Burley AIAA 2001-2210 style.
    H_b = 10^(0.1*(G4(h_over_deltastar_avg, Psi) + g5))
    G_bte = (h*(M^5.5)*L*D)/(r_e^2)*H_b
    SPL_blunt = 10*log10(G_bte)
    return SPL_blunt
end

@concrete struct TEBVSSourceElement{TDirect<:AbstractDirectivity,TUInduction} <: AbstractBroadbandSourceElement{TDirect,TUInduction,NoMachCorrection}
    # Speed of sound, m/s.
    c0
    # Kinematic viscosity, m^2/s
    nu
    # Radial/spanwise length of element, m.
    Δr
    # chord length of element, m.
    chord
    # Trailing edge thickness, m.
    h
    # Solid angle between blade surfaces immediately upstream of the trailing edge, deg.
    Psi
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
function TEBVSSourceElement(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
    return TEBVSSourceElement{BrooksBurleyDirectivity,true}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    TEBVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting trailing edge bluntness-vortex shedding (TEBVS) noise using the BPM/Brooks and Burley method, using position and velocity data expressed in a cylindrical coordinate system.

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
- h: trailing edge thickness (m)
- Psi: solid angle between the blade surfaces immediately upstream of the trailing edge (rad)
- vn: normal velocity of fluid (m/s)
- vr: radial velocity of fluid (m/s)
- vc: circumferential velocity of the fluid (m/s)
- τ: source time (s)
- Δτ: source time duration (s)
- bl: Boundary layer struct, i.e. an AbstractBoundaryLayer.
- twist_about_positive_y: if `true`, apply twist ϕ about positive y axis, negative y axis otherwise
"""
function TEBVSSourceElement{TDirect,TUInduction}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) where {TDirect,TUInduction}
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
    return TEBVSSourceElement{TDirect,TUInduction}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, and include induction in the flow speed normal to span (TUInduction == true).
function TEBVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    return TEBVSSourceElement{BrooksBurleyDirectivity,true}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::TEBVSSourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::TEBVSSourceElement{TDirect,TUInduction}) where {TDirect,TUInduction}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return TEBVSSourceElement{TDirect,TUInduction}(se.c0, se.nu, se.Δr, se.chord, se.h, se.Psi, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.chord_cross_span_to_get_top_uvec)
end

function _teb_vs(freq, h_U, h_over_deltastar_avg, St_3pp, Psi, g4, G_teb_vs_scaler)
    # Equation 71 from the BPM report.
    St_3p = freq*h_U
    St_3prime_over_St_3prime_peak = St_3p/St_3pp
    g5temp = G5(h_over_deltastar_avg, Psi, St_3prime_over_St_3prime_peak)

    # This next check is in the code listing in the BPM report appendix.
    # Need to find G5 for h_over_deltastar_avg = 0.25 for the F4TEMP variable.
    f4temp = G5_Psi14(0.25, St_3prime_over_St_3prime_peak)
    g5 = ifelse(g5temp > f4temp, f4temp, g5temp)

    # Equation 70 from the BPM report.
    # SPL_blunt = 10*log10((h*(M^5.5)*L*D)/(r_e^2)) + G4(h_over_deltastar_avg, Psi) + g5
    # Brooks and Burley AIAA 2001-2210 style.
    H_b = 10^(0.1*(g4 + g5))
    G_bte = G_teb_vs_scaler*H_b
    return G_bte
end

function noise(se::TEBVSSourceElement, obs::AbstractAcousticObserver, t_obs, freqs)
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

    # Equation 73 from the BPM report.
    deltastar_avg = 0.5*(deltastar_p + deltastar_s)

    h_over_deltastar_avg = se.h/deltastar_avg
    h_U = se.h/U
    St_3pp = St_3prime_peak(h_over_deltastar_avg, se.Psi)
    g4 = G4(h_over_deltastar_avg, se.Psi)

    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    G_teb_vs_scaler = (se.h*(M^5.5)*se.Δr*Dh)/(r_er^2)
    G_teb_vs = _teb_vs.(freqs, h_U, h_over_deltastar_avg, St_3pp, se.Psi, g4, G_teb_vs_scaler) .* pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)

    # All done.
    return AcousticMetrics.ProportionalBandSpectrumWithTime(G_teb_vs, freqs_obs, dt, t_obs)
end
