# function H_p(St_p, M, Re_c, Re_deltastar_p)

#     return 10^(A(St_p/St1(M), Re_c)/10) * 10^((K_1(Re_c) - 3)/10) * 10^(DeltaK_1(alphastar, Re_deltastar_p)/10)
# end

struct TBLTEBranches
    K_1::Int
    DeltaK_1::Int
    A_min_a0_s::Int
    A_max_a0_s::Int
    A_min_a_s::Int
    A_max_a_s::Int
    A_s::Int
    A_min_a0_p::Int
    A_max_a0_p::Int
    A_min_a_p::Int
    A_max_a_p::Int
    A_p::Int
    B_min_b0::Int
    B_max_b0::Int
    B_min_b::Int
    B_max_b::Int
    B::Int
    St_2::Int
    K_2::Int
end

function K_1_b(Re_c)
    # Equation (47) in the BPM paper.
    if Re_c < 2.47e5
        branch = 1
        return -4.31*log10(Re_c) + 156.3, branch
    elseif Re_c < 8.0e5
        branch = 2
        return -9.0*log10(Re_c) + 181.6, branch
    else
        branch = 3
        return 128.5*one(Re_c), branch
    end
end

function K_1(Re_c)
    # Equation (47) in the BPM paper.
    ret, branch = K_1_b(Re_c)
    return ret
end

function DeltaK_1_b(alphastar, Re_deltastar_p)
    # Equation (48) in the BPM report.
    T = promote_type(typeof(alphastar), typeof(Re_deltastar_p))
    alphastar_deg = alphastar*180/pi
    if Re_deltastar_p < 5000
        branch = 1
        return alphastar_deg*(1.43*log10(Re_deltastar_p) - 5.29), branch
    else
        branch = 2
        return zero(T), branch
    end
end

function DeltaK_1(alphastar, Re_deltastar_p)
    # Equation (48) in the BPM report.
    ret, branch = DeltaK_1_b(alphastar, Re_deltastar_p)
    return ret
end

function St_1(M)
    # Equation (31) from the BPM report.
    return 0.02*M^(-0.6)
end

function A_min_b(a)
    # Equation (35) from the BPM report.
    if a < 0.204
        branch = 1
        return sqrt(67.552 - 886.788*a^2) - 8.219, branch
    elseif a ≤ 0.244
        branch = 2
        return -32.665*a + 3.981, branch
    else
        branch = 3
        return -142.795*a^3 + 103.656*a^2 - 57.757*a + 6.006, branch
    end
end

function A_min(a)
    # Equation (35) from the BPM report.
    ret, branch  = A_min_b(a)
    return ret
end

function A_max_b(a)
    # Equation (36) from the BPM report.
    if a < 0.13
        branch = 1
        return sqrt(67.552 - 886.788*a^2) - 8.219, branch
    elseif a ≤ 0.321
        branch = 2
        return -15.901*a + 1.098, branch
    else
        branch = 3
        return -4.669*a^3 + 3.491*a^2 - 16.699*a + 1.149, branch
    end
end

function A_max(a)
    # Equation (36) from the BPM report.
    ret, branch = A_max_b(a)
    return ret
end

function A_b(St_over_St_peak, Re_c)
    # Equation (37) from the BPM report.
    a = abs(log10(St_over_St_peak))
    
    # Equation (38) from the BPM report.
    if Re_c < 9.52e4
        A_branch = 1
        a0 = 0.57*one(Re_c)
    elseif Re_c ≤ 8.57e5
        A_branch = 2
        a0 = (-9.57e-13)*(Re_c - 8.57e5)^2 + 1.13
    else
        A_branch = 3
        a0 = 1.13*one(Re_c)
    end

    # Equation (39) from the BPM report.
    A_min_a0, A_min_a0_branch = A_min_b(a0)
    A_max_a0, A_max_a0_branch = A_max_b(a0)
    A_R = (-20 - A_min_a0)/(A_max_a0 - A_min_a0)

    # Equation (40) from the BPM report.
    A_min_a, A_min_a_branch = A_min_b(a)
    A_max_a, A_max_a_branch = A_max_b(a)
    return A_min_a + A_R*(A_max_a - A_min_a), A_min_a0_branch, A_max_a0_branch, A_min_a_branch, A_max_a_branch, A_branch
end

function A(St_over_St_peak, Re_c)
    # Equation (37) from the BPM report.
    # Equation (38) from the BPM report.
    # Equation (39) from the BPM report.
    # Equation (40) from the BPM report.
    ret, A_min_a0_branch, A_max_a0_branch, A_min_a_branch, A_max_a_branch, A_branch = A_b(St_over_St_peak, Re_c)
    return ret
end

function B_min_with_branch(b)
    # Equation (41) from the BPM report.
    if b < 0.13
        branch = 1
        return sqrt(16.888 - 886.788*b^2) - 4.109, branch
    elseif b ≤ 0.145
        branch = 2
        return -83.607*b + 8.138, branch
    else
        branch = 3
        return -817.810*b^3 + 355.201*b^2 - 135.024*b + 10.619, branch
    end
end

function B_min(b)
    # Equation (41) from the BPM report.
    ret, branch = B_min_with_branch(b)
    return ret
end

function B_max_with_branch(b)
    # Equation (42) from the BPM report.
    if b < 0.10
        branch = 1
        return sqrt(16.888 - 886.788*b^2) - 4.109, branch
    elseif b ≤ 0.187
        branch = 2
        return -31.330*b + 1.854, branch
    else
        branch = 3
        return -80.541*b^3 + 44.174*b^2 - 39.381*b + 2.344, branch
    end
end

function B_max(b)
    # Equation (42) from the BPM report.
    ret, branch = B_max_with_branch(b)
    return ret
end

function B_b(St_over_St_peak, Re_c)
    # Equation (43) from the BPM report.
    b = abs(log10(St_over_St_peak))

    # Equation (44) from the BPM report.
    if Re_c < 9.52e4
        B_branch = 1
        b0 = 0.30*one(Re_c)
    elseif Re_c ≤ 8.57e5
        B_branch = 2
        b0 = (-4.48e-13)*(Re_c - 8.57e5)^2 + 0.56
    else
        B_branch = 3
        b0 = 0.56*one(Re_c)
    end

    # Equation (45) from the BPM report.
    B_min_b0, B_min_b0_branch = B_min_with_branch(b0)
    B_max_b0, B_max_b0_branch = B_max_with_branch(b0)
    B_R = (-20 - B_min_b0)/(B_max_b0 - B_min_b0)

    # Equation (46) from the BPM report.
    B_min_b, B_min_b_branch = B_min_with_branch(b)
    B_max_b, B_max_b_branch = B_max_with_branch(b)
    return B_min_b + B_R*(B_max_b - B_min_b), B_min_b0_branch, B_max_b0_branch, B_min_b_branch, B_max_b_branch, B_branch
end

function B(St_over_St_peak, Re_c)
    # Equation (43) from the BPM report.
    # Equation (44) from the BPM report.
    # Equation (45) from the BPM report.
    # Equation (46) from the BPM report.
    ret, B_min_b0_branch, B_max_b0_branch, B_min_b_branch, B_max_b_branch, B_branch = B_b(St_over_St_peak, Re_c)
    return ret
end

function St_2_b(St_1, alphastar)
    # Equation (34) from the BPM report.
    T = promote_type(typeof(St_1), typeof(alphastar))
    alphastar_deg = alphastar*180/pi
    if alphastar_deg < 1.333
        branch = 1
        return St_1*one(T), branch
    elseif alphastar_deg ≤ 12.5
        branch = 2
        return St_1*10.0^(0.0054*(alphastar_deg - 1.333)^2), branch
    else
        branch = 3
        return St_1*4.72*one(T), branch
    end
end

function St_2(St_1, alphastar)
    # Equation (34) from the BPM report.
    ret, branch = St_2_b(St_1, alphastar)
    return ret
end

# function gammas_betas(M)
#     # Equation (50) from the BPM report.
#     gamma_deg = 27.094*M + 3.31
#     gamma0_deg = 23.43*M + 4.651
#     beta_deg = 72.65*M + 10.74
#     beta0_deg = -34.19*M - 13.82
#     return gamma_deg, gamma0_deg, beta_deg, beta0_deg
# end

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

function K_2_b(Re_c, M, alphastar)
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
        branch = 1
        return k_1 - 1000, branch
    # elseif alphastar_deg ≤ gamma0_deg + gamma_deg
    elseif alphastar_deg ≤ gamma0_deg + sqrt(-(gamma_deg/beta_deg)^2*(-12 - beta0_deg)^2 + gamma_deg^2)
        branch = 2
        return k_1 + sqrt(beta_deg^2 - (beta_deg/gamma_deg)^2*(alphastar_deg - gamma0_deg)^2) + beta0_deg, branch
    else
        branch = 3
        return k_1 - 12, branch
    end
end

function K_2(Re_c, M, alphastar)
    # Equation (50) from the BPM report.
    # Equation (49) from the BPM report.
    ret, branch = K_2_b(Re_c, M, alphastar)
    return ret
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
    D = Dbar_h(theta_e, phi_e, M, M_c)
    Re_c = U*chord/nu
    deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord

    St_s = freq*deltastar_s/U
    St_peak = St_1(M)
    A_s = A(St_s/St_peak, Re_c)

    SPL_s = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + A_s + K_1(Re_c) - 3
    return SPL_s
end

function TBL_TE_p(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    D = Dbar_h(theta_e, phi_e, M, M_c)
    Re_c = U*chord/nu
    deltastar_p = disp_thickness_p(bl, Re_c, alphastar)*chord

    St_p = freq*deltastar_p/U
    St_peak = St_1(M)
    A_p = A(St_p/St_peak, Re_c)

    Re_deltastar_p = U*deltastar_p/nu
    ΔK_1 = DeltaK_1(alphastar, Re_deltastar_p)

    SPL_p = 10*log10((deltastar_p*M^5*L*D)/(r_e^2)) + A_p + K_1(Re_c) - 3 + ΔK_1
    return SPL_p
end

function TBL_TE_alpha(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, bl)
    D = Dbar_h(theta_e, phi_e, M, M_c)
    Re_c = U*chord/nu
    deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord

    St_s = freq*deltastar_s/U
    St_peak = St_2(St_1(M), alphastar)

    SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + B(St_s/St_peak, Re_c) + K_2(Re_c, M, alphastar)
    return SPL_alpha
end

function TBL_TE_branch(freq, nu, L, chord, U, M, M_c, r_e, theta_e, phi_e, alphastar, alphastar0, bl)
    T = promote_type(typeof(freq), typeof(nu), typeof(L), typeof(chord), typeof(U), typeof(M), typeof(M_c), typeof(r_e), typeof(theta_e), typeof(phi_e), typeof(alphastar), typeof(alphastar0))
    Re_c = U*chord/nu

    # gamma_deg, gamma0_deg, beta_deg, beta0_deg = gammas_betas(M)
    gamma0_deg = gamma0(M)
    deltastar_s = disp_thickness_s(bl, Re_c, alphastar)*chord
    St_s = freq*deltastar_s/U

    if alphastar*180/pi > min(gamma0_deg, alphastar0*180/pi)
        SPL_s = -100*one(T)
        SPL_p = -100*one(T)
        D = Dbar_l(theta_e, phi_e, M)
        St_peak_p = St_1(M)
        St_peak_alpha, St_2_branch = St_2_b(St_peak_p, alphastar)
        A_prime, A_min_a0_s_branch, A_max_a0_s_branch, A_min_a_s_branch, A_max_a_s_branch, A_s_branch = A_b(St_s/St_peak_alpha, 3*Re_c)
        k2, K_2_branch = K_2_b(Re_c, M, alphastar)
        SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + A_prime + k2

        K_1_branch = 0
        DeltaK_1_branch = 0
        A_min_a0_p_branch = 0
        A_max_a0_p_branch = 0
        A_min_a_p_branch = 0
        A_max_a_p_branch = 0
        A_p_branch = 0
        B_min_b0_branch = 0
        B_max_b0_branch = 0
        B_min_b_branch = 0
        B_max_b_branch = 0
        B_branch = 0
    else
        D = Dbar_h(theta_e, phi_e, M, M_c)

        St_peak_p = St_1(M)
        St_peak_alpha, St_2_branch = St_2_b(St_peak_p, alphastar)
        St_peak_s = 0.5*(St_peak_p + St_peak_alpha)

        A_s, A_min_a0_s_branch, A_max_a0_s_branch, A_min_a_s_branch, A_max_a_s_branch, A_s_branch = A_b(St_s/St_peak_s, Re_c)

        k_1, K_1_branch = K_1_b(Re_c)
        SPL_s = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + A_s + k_1 - 3

        deltastar_p = disp_thickness_p(bl, Re_c, alphastar)*chord

        St_p = freq*deltastar_p/U
        A_p, A_min_a0_p_branch, A_max_a0_p_branch, A_min_a_p_branch, A_max_a_p_branch, A_p_branch = A_b(St_p/St_peak_p, Re_c)

        Re_deltastar_p = U*deltastar_p/nu
        ΔK_1, DeltaK_1_branch = DeltaK_1_b(alphastar, Re_deltastar_p)

        SPL_p = 10*log10((deltastar_p*M^5*L*D)/(r_e^2)) + A_p + k_1 - 3 + ΔK_1

        B_alpha, B_min_b0_branch, B_max_b0_branch, B_min_b_branch, B_max_b_branch, B_branch = B_b(St_s/St_peak_alpha, Re_c)
        k2, K_2_branch = K_2_b(Re_c, M, alphastar)
        SPL_alpha = 10*log10((deltastar_s*M^5*L*D)/(r_e^2)) + B_alpha + k2
    end

    tblte_branches = TBLTEBranches(
        K_1_branch,
        DeltaK_1_branch,
        A_min_a0_s_branch,
        A_max_a0_s_branch,
        A_min_a_s_branch,
        A_max_a_s_branch,
        A_s_branch,
        A_min_a0_p_branch,
        A_max_a0_p_branch,
        A_min_a_p_branch,
        A_max_a_p_branch,
        A_p_branch,
        B_min_b0_branch,
        B_max_b0_branch,
        B_min_b_branch,
        B_max_b_branch,
        B_branch,
        St_2_branch,
        K_2_branch)
    return SPL_s, SPL_p, SPL_alpha, tblte_branches
end
