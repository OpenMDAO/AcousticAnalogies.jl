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

abstract type AbstractMachCorrection end
struct NoMachCorrection <: AbstractMachCorrection end
struct PrandtlGlauertMachCorrection <: AbstractMachCorrection end

struct TBLTESourceElement{
        TDirect<:AbstractDirectivity,TUInduction,TMachCorrection,TDoppler,
        Tc0,Tnu,TΔr,Tchord,Ty0dot,Ty1dot,Ty1dot_fluid,Tτ,TΔτ,Tspan_uvec,Tchord_uvec,Tbl
    } <: AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}
    # Speed of sound, m/s.
    c0::Tc0
    # Kinematic viscosity, m^2/s
    nu::Tnu
    # Radial/spanwise length of element, m.
    Δr::TΔr
    # chord length of element, m.
    chord::Tchord
    # Source position, m.
    y0dot::Ty0dot
    # Source velocity, m/s.
    y1dot::Ty1dot
    # Fluid velocity, m/s.
    y1dot_fluid::Ty1dot_fluid
    # Source time, s.
    τ::Tτ
    # Time step size, i.e. the amount of time this source element "exists" with these properties, s.
    Δτ::TΔτ
    # Radial/spanwise unit vector, aka unit vector aligned with the element's span direction.
    span_uvec::Tspan_uvec
    # Chordwise unit vector, aka unit vector aligned with the element's chord line, pointing from leading edge to trailing edge.
    chord_uvec::Tchord_uvec
    # Boundary layer struct, i.e. an AbstractBoundaryLayer.
    bl::Tbl
    # `Bool` indicating chord_uvec×span_uvec will give a vector pointing from bottom side (usually pressure side) to top side (usually suction side) if `true`, or the opposite if `false`.
    chord_cross_span_to_get_top_uvec::Bool

    function TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, Δr, chord, y0dot::AbstractVector, y1dot::AbstractVector, y1dot_fluid::AbstractVector, τ, Δτ, span_uvec::AbstractVector, chord_uvec::AbstractVector, bl, chord_cross_span_to_get_top_uvec::Bool) where {TDirect<:AbstractDirectivity,TUInduction,TMachCorrection,TDoppler}
        return new{
            TDirect,TUInduction,TMachCorrection,TDoppler,
            typeof(c0), typeof(nu), typeof(Δr), typeof(chord), typeof(y0dot), typeof(y1dot), typeof(y1dot_fluid), typeof(τ), typeof(Δτ), typeof(span_uvec), typeof(chord_uvec), typeof(bl)
        }(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
    end
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the PrandtlGlauertMachCorrection, and Doppler-shift.
function TBLTESourceElement(c0, nu, Δr, chord, y0dot::AbstractVector, y1dot::AbstractVector, y1dot_fluid::AbstractVector, τ, Δτ, span_uvec::AbstractVector, chord_uvec::AbstractVector, bl, chord_cross_span_to_get_top_uvec::Bool)
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
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
function TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y::Bool) where {TDirect,TUInduction,TMachCorrection,TDoppler}
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
    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec::Bool)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the PrandtlGlauertMachCorrection, and Doppler-shift.
function TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y::Bool)
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
end

"""
    TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, U, α, τ, Δτ, bl, twist_about_positive_y)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE) noise using the BPM/Brooks and Burley method, using the velocity magnitude `U` and angle of attack `α`.

The `r` and `θ` arguments are used to define the radial and circumferential position of the source element in a cylindrical coordinate system.
The `U` and `α` arguments are the velocity magnitude normal to the source element length and the angle of attack, respectively.
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
- U: velocity magnitude (m/s)
- α: angle of attack (rad)
- τ: source time (s)
- Δτ: source time duration (s)
- bl: Boundary layer struct, i.e. an AbstractBoundaryLayer.
- twist_about_positive_y: if `true`, apply twist ϕ about positive y axis, negative y axis otherwise
"""
function TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, r, θ, Δr, chord, ϕ, U, α, τ, Δτ, bl, twist_about_positive_y) where {TDirect,TUInduction,TMachCorrection,TDoppler}
    precone = 0
    pitch = 0
    phi = ϕ - α
    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(ϕ, precone, pitch, r, θ, U, phi, twist_about_positive_y)
    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

# Default to using the `BrooksBurleyDirectivity` directivity function, include induction in the flow speed normal to span (TUInduction == true), use the PrandtlGlauertMachCorrection, and Doppler-shift.
function TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, U, α, τ, Δτ, bl, twist_about_positive_y::Bool)
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(c0, nu, r, θ, Δr, chord, ϕ, U, α, τ, Δτ, bl, twist_about_positive_y)
end

"""
    (trans::KinematicTransformation)(se::TBLTESourceElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}) where {TDirect,TUInduction,TMachCorrection,TDoppler}
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)
    y0dot, y1dot_fluid = trans(se.τ, se.y0dot, se.y1dot_fluid, linear_only)
    linear_only = true
    span_uvec = trans(se.τ, se.span_uvec, linear_only)
    chord_uvec = trans(se.τ, se.chord_uvec, linear_only)

    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(se.c0, se.nu, se.Δr, se.chord, y0dot, y1dot, y1dot_fluid, se.τ, se.Δτ, span_uvec, chord_uvec, se.bl, se.chord_cross_span_to_get_top_uvec)
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
