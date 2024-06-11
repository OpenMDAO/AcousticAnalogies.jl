abstract type AbstractTipAlphaCorrection end

struct NoTipAlphaCorrection <: AbstractTipAlphaCorrection end
struct BPMTipAlphaCorrection <: AbstractTipAlphaCorrection end

abstract type AbstractBladeTip{TTipAlphaCorrection} end

alpha_zerolift(blade_tip::AbstractBladeTip) = blade_tip.alpha0lift
tip_alpha_correction(blade_tip::AbstractBladeTip) = blade_tip.tip_alpha_correction

struct RoundedTip{TTipAlphaCorrection,TAlpha0Lift} <: AbstractBladeTip{TTipAlphaCorrection}
    tip_alpha_correction::TTipAlphaCorrection
    alpha0lift::TAlpha0Lift

    function RoundedTip(tip_alpha_correction::TTipAlphaCorrection, alpha0lift::TAlpha0Lift) where {TTipAlphaCorrection<:AbstractTipAlphaCorrection,TAlpha0Lift}
        return new{TTipAlphaCorrection, TAlpha0Lift}(tip_alpha_correction, alpha0lift)
    end
end
RoundedTip(alpha0lift=0.0) = RoundedTip(NoTipAlphaCorrection(), alpha0lift)

struct FlatTip{TTipAlphaCorrection,TAlpha0Lift} <: AbstractBladeTip{TTipAlphaCorrection}
    tip_alpha_correction::TTipAlphaCorrection
    alpha0lift::TAlpha0Lift

    function FlatTip(tip_alpha_correction::TTipAlphaCorrection, alpha0lift::TAlpha0Lift) where {TTipAlphaCorrection<:AbstractTipAlphaCorrection,TAlpha0Lift}
        return new{TTipAlphaCorrection, TAlpha0Lift}(tip_alpha_correction, alpha0lift)
    end
end
FlatTip(alpha0lift=0.0) = FlatTip(NoTipAlphaCorrection(), alpha0lift)

function tip_vortex_alpha_correction(blade_tip::AbstractBladeTip{NoTipAlphaCorrection}, alphatip)
    return alphatip
end

function tip_vortex_alpha_correction(blade_tip::AbstractBladeTip{BPMTipAlphaCorrection}, alphatip)
    # Referencing the tip vortex angle of attack correction to the zero-lift angle of attack ensures that the correction will never cause the angle of attack to drop below the zero-lift value, assuming the correction factor (0.71 here) is between 0 and 1.
    return 0.71*(alphatip - alpha_zerolift(blade_tip)) + alpha_zerolift(blade_tip)
    # return 0.71*alphatip + (1 - 0.71)*alpha_zerolift(blade_tip)
end

function tip_vortex_size_c(::RoundedTip, alphatip)
    # Equation 63 in the BPM report.
    alphatip_deg = alphatip*180/pi
    return 0.008*alphatip_deg
end

function tip_vortex_size_c(::FlatTip, alphatip)
    alphatip_deg = alphatip*180/pi
    # Equation 67 in the BPM report.
    if alphatip_deg < 2
        return 0.0230 + 0.0169*alphatip_deg
    else
        return 0.0378 + 0.0095*alphatip_deg
    end
end

function tip_vortex_max_mach_number(::AbstractBladeTip, M, alphatip)
    alphatip_deg = alphatip*180/pi
    # Equation 64 in the BPM report.
    M_max = (1 + 0.036*alphatip_deg)*M
    return M_max
end

function TIP(freq, chord, M, M_c, U_max, M_max, r_e, theta_e, phi_e, alphatip, blade_tip::AbstractBladeTip)
    l = tip_vortex_size_c(blade_tip, alphatip) * chord
    # Equation 62 in the BPM report.
    St_pp = freq*l/U_max

    D = Dbar_h(theta_e, phi_e, M, M_c)

    # Equation 61 in the BPM report.
    # SPL_tip = 10*log10(M^2*M_max^3*l^2*D/r_e^2) - 30.5*(log10(St_pp) + 0.3)^2 + 126
    # Brooks and Burley AIAA 2001-2210 style.
    H_t = 10^(0.1*(-30.5*(log10(St_pp) + 0.3)^2 + 126))
    G_tip = (M^2*M_max^3*l^2*D/r_e^2)*H_t
    SPL_tip = 10*log10(G_tip)
    return SPL_tip
end

@concrete struct TipVortexSourceElement{TDirect<:AbstractDirectivity,TUInduction} <: AbstractBroadbandSourceElement{TDirect,TUInduction,NoMachCorrection}
    # Speed of sound, m/s.
    c0
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
    # Blade tip struct, i.e. an AbstractBladeTip.
    blade_tip
    # `Bool` indicating chord_uvec×span_uvec will give a vector pointing from bottom side (usually pressure side) to top side (usually suction side) if `true`, or the opposite if `false`.
    chord_cross_span_to_get_top_uvec
end

# Default to using the `BrooksBurleyDirectivity` directivity function, and include induction in the flow speed normal to span (TUInduction == true).
function TipVortexSourceElement(c0, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
    return TipVortexSourceElement{BrooksBurleyDirectivity,true}(c0, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

"""
    TipVortexSourceElement(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)

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
- blade_tip: Blade tip struct, i.e. an AbstractBladeTip.
- twist_about_positive_y: if `true`, apply twist ϕ about positive y axis, negative y axis otherwise
"""
function TipVortexSourceElement{TDirect,TUInduction}(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) where {TDirect,TUInduction}
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
    return TipVortexSourceElement{TDirect,TUInduction}(c0, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, and include induction in the flow speed normal to span (TUInduction == true).
function TipVortexSourceElement(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
    return TipVortexSourceElement{BrooksBurleyDirectivity,true}(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::TipVortexSourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::TipVortexSourceElement{TDirect,TUInduction}) where {TDirect,TUInduction}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return TipVortexSourceElement{TDirect,TUInduction}(se.c0, se.Δr, se.chord, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.blade_tip, se.chord_cross_span_to_get_top_uvec)
end

function _tip(freq, l_U_max, scaler)
    St_pp = freq*l_U_max
    H_t = 10^(0.1*(-30.5*(log10(St_pp) + 0.3)^2 + 126))
    # G_tip = (M^2*M_max^3*l^2*D/r_e^2)*H_t
    G_tip = scaler*H_t
    return G_tip
end

function noise(se::TipVortexSourceElement, obs::AbstractAcousticObserver, t_obs, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
    # Position of the observer:
    x_obs = obs(t_obs)

    # Need the angle of attack, including the possible tip correction.
    alphatip = tip_vortex_alpha_correction(se.blade_tip, angle_of_attack(se))

    # Need the directivity functions.
    top_is_suction = is_top_suction(se.bl, alphatip)
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

    # Now that we've decided on the directivity functions and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphatip_positive = abs_cs_safe(alphatip - alpha_zerolift(se.bl))

    # Mach number of the flow speed normal to span.
    M = U/se.c0

    # Need the maximum mach number near the tip vortex.
    M_max = tip_vortex_max_mach_number(se.blade_tip, M, alphatip_positive)

    # Now we can find the maximum speed near the tip vortex.
    U_max = M_max * se.c0

    # Get the tip vortex size.
    l = tip_vortex_size_c(se.blade_tip, alphatip_positive) * se.chord

    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    l_U_max = l/U_max
    G_tip_scaler = (M^2*M_max^3*l^2*Dh/r_er^2)
    G_tip = _tip.(freqs, l_U_max, G_tip_scaler) .* pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)

    # All done.
    return AcousticMetrics.ProportionalBandSpectrumWithTime(G_tip, freqs_obs, dt, t_obs)
end
