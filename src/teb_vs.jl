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
    deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord
    deltastar_p = disp_thickness_p(bl, Re_c, alphastar)*chord

    # Equation 73 from the BPM report.
    deltastar_avg = 0.5*(deltastar_p + deltastar_s)

    h_over_deltastar_avg = h/deltastar_avg
    D = Dbar_h(theta_e, phi_e, M, M_c)

    # Equation 71 from the BPM report.
    St_3p = freq*h/U

    St_3prime_over_St_3prime_peak = St_3p/St_3prime_peak(h_over_deltastar_avg, Psi)

    g5 = G5(h_over_deltastar_avg, Psi, St_3prime_over_St_3prime_peak)

    # This next check is in the code listing in the BPM report appendix.
    # Need to find G5 for h_over_deltastar_avg = 0.25 for the F4TEMP variable.
    f4temp = G5_Psi14(0.25, St_3prime_over_St_3prime_peak)
    if g5 > f4temp
        g5 = f4temp
    end

    # Equation 70 from the BPM report.
    # SPL_blunt = 10*log10((h*(M^5.5)*L*D)/(r_e^2)) + G4(h_over_deltastar_avg, Psi) + g5
    # Brooks and Burley AIAA 2001-2210 style.
    H_b = 10^(0.1*(G4(h_over_deltastar_avg, Psi) + g5))
    G_bte = (h*(M^5.5)*L*D)/(r_e^2)*H_b
    SPL_blunt = 10*log10(G_bte)
    return SPL_blunt
end
