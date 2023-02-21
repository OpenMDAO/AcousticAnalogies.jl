function H_p(St_p, M, Re_c, Re_deltastar_p)

    return 10^(A(St_p/St1(M), Re_c)/10) * 10^((K_1(Re_c) - 3)/10) * 10^(DeltaK_1(alphastar, Re_deltastar_p)/10)
end

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
    if Re_deltastar_p < 50000
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
    T = promote_type(typeof(St_1), typeof(alphastar))
    alphastar_deg = alphastar*180/pi
    if alphastar_deg < 1.33
        return St_1*one(T)
    elseif alphastar_deg ≤ 12.5
        return St_1*10.0^(0.0054*(alphastar_deg - 1.33)^2)
    else
        return St_1*4.72*one(T)
    end
end

function K_2(Re_c, M, alphastar)
    T = promote_type(typeof(Re_c), typeof(M), typeof(alphastar))
    alphastar_deg = alphastar*180/pi

    k_1 = K_1(Re_c)*one(T)
    gamma = 27.094*M + 3.31
    gamma0 = 23.43*M + 4.651
    beta = 72.65*M + 10.74
    beta0 = -34.19*M - 13.82

    if alphastar_deg < gamma0 - gamma
        return k_1 - 1000
    # elseif alphastar_deg ≤ gamma0 + gamma
    elseif alphastar_deg ≤ gamma0 + sqrt(-(gamma/beta)^2*(-12 - beta0)^2 + gamma^2)
        return k_1 + sqrt(beta^2 - (beta/gamma)^2*(alphastar_deg - gamma0)^2) + beta0
    else
        return k_1 - 12
    end
end
