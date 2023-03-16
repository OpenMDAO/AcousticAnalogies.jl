function St_3prime_peak(h_over_deltastar_avg, Phi)
    Phi_deg = Phi*180/pi
    # Equation 72 from the BPM report.
    if h_over_deltastar_avg < 0.2
        return 0.1*(h_over_deltastar_avg) + 0.095 - 0.00243*Phi_deg
    else
        return (0.212 - 0.0045*Phi_deg)/(1 + 0.235/h_over_deltastar_avg - 0.0132/(h_over_deltastar_avg)^2)
    end
end

function G4(h_over_deltastar_avg, Phi)
    Phi_deg = Phi*180/pi
    # Equation 74 from the BPM report.
    if h_over_deltastar_avg ≤ 5
        return 17.5*log10(h_over_deltastar_avg) + 157.5 - 1.114*Phi_deg
    else
        return 169.7 - 1.114*Phi_deg
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
