@concrete struct CombinedNoTipBroadbandSourceElement{TDirect<:AbstractDirectivity,TUInduction,TMachCorrection,TDoppler} <: AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}
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
    # Solid angle between blade surfaces immediately upstream of the trailing edge, rad.
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

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the Prandtl-Glauert mach number correction, and Doppler-shift.
function CombinedNoTipBroadbandSourceElement(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
    return CombinedNoTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    CombinedNoTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE), laminar boundary layer-vortex shedding (LBLVS) noise, and trailing edge bluntness-vortex shedding (TEBVS) noise using the BPM/Brooks and Burley method, using position and velocity data expressed in a cylindrical coordinate system.

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
function CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) where {TDirect,TUInduction,TMachCorrection,TDoppler}
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
    return CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the Prandtl-Glauert mach number correction, and Doppler-shift.
function CombinedNoTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
    return CombinedNoTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::CombinedNoTipBroadbandSourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}) where {TDirect,TUInduction,TMachCorrection,TDoppler}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(se.c0, se.nu, se.Δr, se.chord, se.h, se.Psi, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.chord_cross_span_to_get_top_uvec)
end

"""
    CombinedNoTipOutput(G_s, G_p, G_alpha, G_teb, cbands, dt, t)

Output of the combined broadband noise calculation not including tip vortex noise: the acoustic pressure autospectrum centered at time `t` over observer duration `dt` and observer frequencies `cbands` for the TBLTE suction side `G_s`, TBLTE pressure side `G_p`, TBLTE separation noise `G_alpha`, and trailing edge bluntness noise `G_teb`.
"""
struct CombinedNoTipOutput{NO,TF,TG<:AbstractVector{TF},TFreqs<:AcousticMetrics.AbstractProportionalBands{NO,:center},TDTime,TTime} <: AcousticMetrics.AbstractProportionalBandSpectrum{NO,TF}
    G_s::TG
    G_p::TG
    G_alpha::TG
    G_lblvs::TG
    G_teb::TG
    cbands::TFreqs
    dt::TDTime
    t::TTime

    function CombinedNoTipOutput(G_s::TG, G_p::TG, G_alpha::TG, G_lblvs, G_teb::TG, cbands::AcousticMetrics.AbstractProportionalBands{NO,:center}, dt, t) where {NO,TG}
        ncbands = length(cbands)
        length(G_s) == ncbands || throw(ArgumentError("length(G_s) must match length(cbands)"))
        length(G_p) == ncbands || throw(ArgumentError("length(G_p) must match length(cbands)"))
        length(G_alpha) == ncbands || throw(ArgumentError("length(G_alpha) must match length(cbands)"))
        length(G_lblvs) == ncbands || throw(ArgumentError("length(G_lblvs) must match length(cbands)"))
        length(G_teb) == ncbands || throw(ArgumentError("length(G_teb) must match length(cbands)"))
        dt > zero(dt) || throw(ArgumentError("dt must be positive"))
        return new{NO,eltype(TG),TG,typeof(cbands),typeof(dt),typeof(t)}(G_s, G_p, G_alpha, G_lblvs, G_teb, cbands, dt, t)
    end
end

@inline function Base.getindex(pbs::CombinedNoTipOutput, i::Int)
    @boundscheck checkbounds(pbs, i)
    return @inbounds pbs.G_s[i] + pbs.G_p[i] + pbs.G_alpha[i] + +pbs.G_lblvs[i] + pbs.G_teb[i]
end

@inline AcousticMetrics.has_observer_time(pbs::CombinedNoTipOutput) = true
@inline AcousticMetrics.observer_time(pbs::CombinedNoTipOutput) = pbs.t
@inline AcousticMetrics.timestep(pbs::CombinedNoTipOutput) = pbs.dt
@inline AcousticMetrics.time_scaler(pbs::CombinedNoTipOutput, period) = timestep(pbs)/period

function noise(se::CombinedNoTipBroadbandSourceElement, obs::AbstractAcousticObserver, t_obs, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
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

    # Need the boundary layer thickness for the pressure side for LBL-VS noise.
    delta_p = bl_thickness_p(se.bl, Re_c, alphastar)*se.chord

    # Now that we've decided on the directivity functions and the displacement thickness, and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphastar_positive = abs_cs_safe(alphastar - alpha_zerolift(se.bl))

    # Mach number of the flow speed normal to span.
    M = U/se.c0

    # This stuff is used to decide if the blade element is stalled or not.
    alphastar0 = alpha_stall(se.bl, Re_c)
    gamma0_deg = gamma0(M)
    deep_stall = (alphastar_positive*180/pi) > min(gamma0_deg, alphastar0*180/pi)
    # if deep_stall
    #     println("deep_stall! M = $(M), alphastar_positive*180/pi = $(alphastar_positive*180/pi), gamma0_deg = $(gamma0_deg), alphastar0*180/pi = $(alphastar0*180/pi)")
    #     println("forcing deep_stall == false")
    #     deep_stall = false
    # end

    St_peak_p = St_1(M)
    St_peak_alpha = St_2(St_peak_p, alphastar_positive)
    St_peak_s = 0.5*(St_peak_p + St_peak_alpha)

    Re_deltastar_p = U*deltastar_p/se.nu
    k_1 = K_1(Re_c)
    k_2 = K_2(Re_c, M, alphastar_positive)
    Δk_1 = DeltaK_1(alphastar_positive, Re_deltastar_p)

    deltastar_s_U = deltastar_s/U
    deltastar_p_U = deltastar_p/U

    # Stuff for LBLVS noise.
    delta_p_U = delta_p/U
    St_p_p = St_peak_prime(St_1_prime(Re_c), alphastar_positive)
    Re_c_over_Re_c0 = Re_c / Re_c0(alphastar_positive)
    g2 = G2(Re_c_over_Re_c0)
    g3 = G3(alphastar_positive)

    # Brooks and Burley 2001 recommend a Prandtl-Glauert style Mach number correction, but only for the TBLTE noise.
    # But whether or not it's included is dependent on the TMachCorrection type parameter for the source element.
    m_corr = mach_correction(se, M)

    # Equation 73 from the BPM report.
    deltastar_avg = 0.5*(deltastar_p + deltastar_s)

    h_over_deltastar_avg = se.h/deltastar_avg
    h_U = se.h/U
    St_3pp = St_3prime_peak(h_over_deltastar_avg, se.Psi)
    g4 = G4(h_over_deltastar_avg, se.Psi)

    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    G_s_scaler = (deltastar_s*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_s = _tble_te_s.(freqs, deltastar_s_U, Re_c, St_peak_s, k_1, G_s_scaler, deep_stall).*pref2

    G_p_scaler = (deltastar_p*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_p = _tble_te_p.(freqs, deltastar_p_U, Re_c, St_peak_p, k_1, Δk_1, G_p_scaler, deep_stall).*pref2

    G_alpha_scaler_l = (deltastar_s*M^5*se.Δr*Dl)/(r_er^2)*m_corr
    G_alpha_scaler_h = G_s_scaler
    G_alpha = _tble_te_alpha.(freqs, Re_c, deltastar_s_U, St_peak_alpha, k_2, G_alpha_scaler_l, G_alpha_scaler_h, deep_stall).*pref2

    G_lbl_vs_scaler = (delta_p*M^5*se.Δr*Dh)/(r_er^2)
    G_lbl_vs = _lbl_vs.(freqs, delta_p_U, St_p_p, g2, g3, G_lbl_vs_scaler) .* pref2

    G_teb_vs_scaler = (se.h*(M^5.5)*se.Δr*Dh)/(r_er^2)
    G_teb_vs = _teb_vs.(freqs, h_U, h_over_deltastar_avg, St_3pp, se.Psi, g4, G_teb_vs_scaler) .* pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)
    @assert AcousticMetrics.freq_scaler(freqs_obs) ≈ doppler

    # All done.
    return CombinedNoTipOutput(G_s, G_p, G_alpha, G_lbl_vs, G_teb_vs, freqs_obs, dt, t_obs)
end

@concrete struct CombinedWithTipBroadbandSourceElement{TDirect<:AbstractDirectivity,TUInduction,TMachCorrection,TDoppler} <: AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}
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
    # Solid angle between blade surfaces immediately upstream of the trailing edge, rad.
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
    # Blade tip struct, i.e. an AbstractBladeTip.
    blade_tip
    # `Bool` indicating chord_uvec×span_uvec will give a vector pointing from bottom side (usually pressure side) to top side (usually suction side) if `true`, or the opposite if `false`.
    chord_cross_span_to_get_top_uvec
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the Prandtl-Glauert mach number correction, and Doppler-shift.
function CombinedWithTipBroadbandSourceElement(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
    return CombinedWithTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

"""
    CombinedWithTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE), laminar boundary layer-vortex shedding (LBLVS) noise, trailing edge bluntness-vortex shedding (TEBVS) noise, and tip vortex noise using the BPM/Brooks and Burley method, using position and velocity data expressed in a cylindrical coordinate system.

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
- blade_tip: Blade tip struct, i.e. an AbstractBladeTip.
- twist_about_positive_y: if `true`, apply twist ϕ about positive y axis, negative y axis otherwise
"""
function CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) where {TDirect,TUInduction,TMachCorrection,TDoppler}
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
    return CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, Δr, chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the Prandtl-Glauert mach number correction, and Doppler-shift.
function CombinedWithTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
    return CombinedWithTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::CombinedWithTipBroadbandSourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}) where {TDirect,TUInduction,TMachCorrection,TDoppler}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(se.c0, se.nu, se.Δr, se.chord, se.h, se.Psi, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.blade_tip, se.chord_cross_span_to_get_top_uvec)
end

"""
    CombinedWithTipOutput(G_s, G_p, G_alpha, G_teb, G_tip, cbands, dt, t)

Output of the combined broadband noise calculation: the acoustic pressure autospectrum centered at time `t` over observer duration `dt` and observer frequencies `cbands` for the TBLTE suction side `G_s`, TBLTE pressure side `G_p`, TBLTE separation noise `G_alpha`, trailing edge bluntness noise `G_teb`, and tip vortex noise `G_tip`.
"""
struct CombinedWithTipOutput{NO,TF,TG<:AbstractVector{TF},TFreqs<:AcousticMetrics.AbstractProportionalBands{NO,:center},TDTime,TTime} <: AcousticMetrics.AbstractProportionalBandSpectrum{NO,TF}
    G_s::TG
    G_p::TG
    G_alpha::TG
    G_lblvs::TG
    G_teb::TG
    G_tip::TG
    cbands::TFreqs
    dt::TDTime
    t::TTime

    function CombinedWithTipOutput(G_s::TG, G_p::TG, G_alpha::TG, G_lblvs, G_teb::TG, G_tip::TG, cbands::AcousticMetrics.AbstractProportionalBands{NO,:center}, dt, t) where {NO,TG}
        ncbands = length(cbands)
        length(G_s) == ncbands || throw(ArgumentError("length(G_s) must match length(cbands)"))
        length(G_p) == ncbands || throw(ArgumentError("length(G_p) must match length(cbands)"))
        length(G_alpha) == ncbands || throw(ArgumentError("length(G_alpha) must match length(cbands)"))
        length(G_lblvs) == ncbands || throw(ArgumentError("length(G_lblvs) must match length(cbands)"))
        length(G_teb) == ncbands || throw(ArgumentError("length(G_teb) must match length(cbands)"))
        length(G_tip) == ncbands || throw(ArgumentError("length(G_tip) must match length(cbands)"))
        dt > zero(dt) || throw(ArgumentError("dt must be positive"))
        return new{NO,eltype(TG),TG,typeof(cbands),typeof(dt),typeof(t)}(G_s, G_p, G_alpha, G_lblvs, G_teb, G_tip, cbands, dt, t)
    end
end

@inline function Base.getindex(pbs::CombinedWithTipOutput, i::Int)
    @boundscheck checkbounds(pbs, i)
    return @inbounds pbs.G_s[i] + pbs.G_p[i] + pbs.G_alpha[i] + pbs.G_teb[i] + pbs.G_tip[i]
end

@inline AcousticMetrics.has_observer_time(pbs::CombinedWithTipOutput) = true
@inline AcousticMetrics.observer_time(pbs::CombinedWithTipOutput) = pbs.t
@inline AcousticMetrics.timestep(pbs::CombinedWithTipOutput) = pbs.dt
@inline AcousticMetrics.time_scaler(pbs::CombinedWithTipOutput, period) = timestep(pbs)/period

function noise(se::CombinedWithTipBroadbandSourceElement, obs::AbstractAcousticObserver, t_obs, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
    # Position of the observer:
    x_obs = obs(t_obs)

    # Need the angle of attack.
    alphastar = angle_of_attack(se)

    # Need the angle of attack, including the possible tip correction.
    alphatip = tip_vortex_alpha_correction(se.blade_tip, alphastar)

    # Need the directivity functions.
    # But we have two angles of attack: one that includes a tip vortex correction, one that doesn't.
    # But the angle of attack is only used to determine the value of `top_is_suction`, which in turn only depends on if the angle of attack is greater or less than the zero-lift angle of attack.
    # And the tip vortex alpha correction is designed to not change that.
    # So no worries about which alpha to use.
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

    # Need the boundary layer thickness for the pressure side for LBL-VS noise.
    delta_p = bl_thickness_p(se.bl, Re_c, alphastar)*se.chord

    # Now that we've decided on the directivity functions and the displacement thickness, and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphastar_positive = abs_cs_safe(alphastar - alpha_zerolift(se.bl))

    # Now that we've decided on the directivity functions and we know the correct value of `top_is_suction` we should be able to switch the sign on `alphastar` if it's negative, and reference it to the zero-lift value, as the BPM report does.
    alphatip_positive = abs_cs_safe(alphatip - alpha_zerolift(se.bl))

    # Mach number of the flow speed normal to span.
    M = U/se.c0

    # This stuff is used to decide if the blade element is stalled or not.
    alphastar0 = alpha_stall(se.bl, Re_c)
    gamma0_deg = gamma0(M)
    deep_stall = (alphastar_positive*180/pi) > min(gamma0_deg, alphastar0*180/pi)
    # if deep_stall
    #     println("deep_stall! M = $(M), alphastar_positive*180/pi = $(alphastar_positive*180/pi), gamma0_deg = $(gamma0_deg), alphastar0*180/pi = $(alphastar0*180/pi)")
    #     println("forcing deep_stall == false")
    #     deep_stall = false
    # end

    St_peak_p = St_1(M)
    St_peak_alpha = St_2(St_peak_p, alphastar_positive)
    St_peak_s = 0.5*(St_peak_p + St_peak_alpha)

    Re_deltastar_p = U*deltastar_p/se.nu
    k_1 = K_1(Re_c)
    k_2 = K_2(Re_c, M, alphastar_positive)
    Δk_1 = DeltaK_1(alphastar_positive, Re_deltastar_p)

    deltastar_s_U = deltastar_s/U
    deltastar_p_U = deltastar_p/U

    # Stuff for LBLVS noise.
    delta_p_U = delta_p/U
    St_p_p = St_peak_prime(St_1_prime(Re_c), alphastar_positive)
    Re_c_over_Re_c0 = Re_c / Re_c0(alphastar_positive)
    g2 = G2(Re_c_over_Re_c0)
    g3 = G3(alphastar_positive)

    # Brooks and Burley 2001 recommend a Prandtl-Glauert style Mach number correction, but only for the TBLTE noise.
    # But whether or not it's included is dependent on the TMachCorrection type parameter for the source element.
    m_corr = mach_correction(se, M)

    # Equation 73 from the BPM report.
    deltastar_avg = 0.5*(deltastar_p + deltastar_s)

    h_over_deltastar_avg = se.h/deltastar_avg
    h_U = se.h/U
    St_3pp = St_3prime_peak(h_over_deltastar_avg, se.Psi)
    g4 = G4(h_over_deltastar_avg, se.Psi)

    # Need the maximum mach number near the tip vortex.
    M_max = tip_vortex_max_mach_number(se.blade_tip, M, alphatip_positive)

    # Now we can find the maximum speed near the tip vortex.
    U_max = M_max * se.c0

    # Get the tip vortex size.
    l = tip_vortex_size_c(se.blade_tip, alphatip_positive) * se.chord

    # The Brooks and Burley autospectrums appear to be scaled by the usual squared reference pressure (20 μPa)^2, but I'd like things in dimensional units, so multiply through by that.
    pref2 = 4e-10
    G_s_scaler = (deltastar_s*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_s = _tble_te_s.(freqs, deltastar_s_U, Re_c, St_peak_s, k_1, G_s_scaler, deep_stall).*pref2

    G_p_scaler = (deltastar_p*M^5*se.Δr*Dh)/(r_er^2)*m_corr
    G_p = _tble_te_p.(freqs, deltastar_p_U, Re_c, St_peak_p, k_1, Δk_1, G_p_scaler, deep_stall).*pref2

    G_alpha_scaler_l = (deltastar_s*M^5*se.Δr*Dl)/(r_er^2)*m_corr
    G_alpha_scaler_h = G_s_scaler
    G_alpha = _tble_te_alpha.(freqs, Re_c, deltastar_s_U, St_peak_alpha, k_2, G_alpha_scaler_l, G_alpha_scaler_h, deep_stall).*pref2

    G_lbl_vs_scaler = (delta_p*M^5*se.Δr*Dh)/(r_er^2)
    G_lbl_vs = _lbl_vs.(freqs, delta_p_U, St_p_p, g2, g3, G_lbl_vs_scaler) .* pref2

    G_teb_vs_scaler = (se.h*(M^5.5)*se.Δr*Dh)/(r_er^2)
    G_teb_vs = _teb_vs.(freqs, h_U, h_over_deltastar_avg, St_3pp, se.Psi, g4, G_teb_vs_scaler) .* pref2

    l_U_max = l/U_max
    G_tip_scaler = (M^2*M_max^3*l^2*Dh/r_er^2)
    G_tip = _tip.(freqs, l_U_max, G_tip_scaler) .* pref2

    # Also need the Doppler shift for this source-observer combination.
    doppler = doppler_factor(se, obs, t_obs)

    # Get the doppler-shifted time step and proportional bands.
    dt = se.Δτ/doppler
    freqs_obs = AcousticMetrics.center_bands(freqs, doppler)
    # @assert AcousticMetrics.freq_scaler(freqs_obs) ≈ doppler

    # All done.
    return CombinedWithTipOutput(G_s, G_p, G_alpha, G_lbl_vs, G_teb_vs, G_tip, freqs_obs, dt, t_obs)
end

function pbs_suction(pbs::Union{TBLTEOutput,CombinedNoTipOutput,CombinedWithTipOutput})
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_s, cbands, dt, t)
end

function pbs_pressure(pbs::Union{TBLTEOutput,CombinedNoTipOutput,CombinedWithTipOutput})
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_p, cbands, dt, t)
end

function pbs_alpha(pbs::Union{TBLTEOutput,CombinedNoTipOutput,CombinedWithTipOutput})
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_alpha, cbands, dt, t)
end

function pbs_lblvs(pbs::Union{TBLTEOutput,CombinedNoTipOutput,CombinedWithTipOutput})
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_lblvs, cbands, dt, t)
end

function pbs_teb(pbs::Union{CombinedNoTipOutput,CombinedWithTipOutput})
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_teb, cbands, dt, t)
end

function pbs_tip(pbs::CombinedWithTipOutput)
    t = AcousticMetrics.observer_time(pbs)
    dt = AcousticMetrics.timestep(pbs)
    cbands = AcousticMetrics.center_bands(pbs)
    return AcousticMetrics.ProportionalBandSpectrumWithTime(pbs.G_tip, cbands, dt, t)
end
