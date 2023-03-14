abstract type AbstractBladeTip end

struct RoundedTip <: AbstractBladeTip end
struct FlatTip <: AbstractBladeTip end

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

function TIP(freq, chord, M, M_c, U_max, M_max, r_e, theta_e, phi_e, alphatip, blade_tip::AbstractBladeTip)
    l = tip_vortex_size_c(blade_tip, alphatip) * chord
    # Equation 62 in the BPM report.
    St_pp = freq*l/U_max

    D = Dbar_h(theta_e, phi_e, M, M_c)

    # Equation 61 in the BPM report.
    SPL_tip = 10*log10(M^2*M_max^3*l^2*D/r_e^2) - 30.5*(log10(St_pp) + 0.3)^2 + 126
    return SPL_tip
end
