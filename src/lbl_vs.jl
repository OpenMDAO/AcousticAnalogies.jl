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
    SPL = 10*log10(delta_p*M^5*L*Dbar_h(theta_e, phi_e, M, M_c)/r_e^2) + G1(St_prime_over_St_peak_prime) + G2(Re_c_over_Re_c0) + G3(alphastar)
    return SPL
end
