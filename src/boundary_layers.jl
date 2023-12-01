abstract type AbstractBoundaryLayer end

@concrete struct TrippedN0012BoundaryLayer <: AbstractBoundaryLayer
    alphastar0
end

TrippedN0012BoundaryLayer() = TrippedN0012BoundaryLayer(12.5*pi/180)

@concrete struct UntrippedN0012BoundaryLayer <: AbstractBoundaryLayer
    alphastar0
end

UntrippedN0012BoundaryLayer() = UntrippedN0012BoundaryLayer(12.5*pi/180)

alpha_stall(bl::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}, Re_c) = bl.alphastar0

alpha_zerolift(bl::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}) = zero(bl.alphastar0)

function bl_thickness_0(::TrippedN0012BoundaryLayer, Re_c)
    # Equation 2 from the BPM report.
    logRe_c = log10(Re_c)
    return 10^(1.892 - 0.9045*logRe_c + 0.0596*logRe_c^2)
end

function disp_thickness_0(::TrippedN0012BoundaryLayer, Re_c)
    # Equation 3 from the BPM report.
    if Re_c ≤ 0.3e6
        return 0.0601*Re_c^(-0.114)
    else
        logRe_c = log10(Re_c)
        return 10^(3.411 - 1.5397*logRe_c + 0.1059*logRe_c^2)
    end
end

function bl_thickness_0(::UntrippedN0012BoundaryLayer, Re_c)
    # Equation 5 from the BPM report.
    logRe_c = log10(Re_c)
    return 10^(1.6569 - 0.9045*logRe_c + 0.0596*logRe_c^2)
end

function disp_thickness_0(::UntrippedN0012BoundaryLayer, Re_c)
    # Equation 6 from the BPM report.
    logRe_c = log10(Re_c)
    return 10^(3.0187 - 1.5397*logRe_c + 0.1059*logRe_c^2)
end

function bl_thickness_p(::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}, alphastar)
    # Equation 8 from the BPM report.
    alphastar_deg = alphastar*180/pi
    return 10^(-0.04175*alphastar_deg + 0.00106*alphastar_deg^2)
end

function disp_thickness_p(::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}, alphastar)
    # Equation 9 from the BPM report.
    alphastar_deg = alphastar*180/pi
    return 10^(-0.0432*alphastar_deg + 0.00113*alphastar_deg^2)
end

# function bl_thickness_s(::TrippedN0012BoundaryLayer, alphastar)
#     # Equation 11 from the BPM report.
#     # The report defines the suction-side boundary layer parameters for alphastar values from 0° to 25°.
#     # But what about negative angles of attack?
#     # The NACA0012 airfoil is symmetric, but if the angle of attack goes negative, I guess the pressure and suction sides would switch.
#     # So I'll check that the angle of attack is always positive here.
#     alphastar_deg = alphastar*180/pi
#     if alphastar_deg < 0
#         throw(DomainError(alphastar, "negative alphastar argument invalid"))
#     elseif alphastar_deg ≤ 5
#         return 10^(0.0311*alphastar)
#     elseif alphastar_deg ≤ 12.5
#         return 0.3468*10^(0.1231*alphastar)
#     elseif alphastar_deg ≤ 25
#         return 5.718*10^(0.0258*alphastar)
#     else
#         # What should I do for angles of attack greater than 25°?
#         # Maybe just keep the same thickness?
#         return 5.718*10^(0.0258*25*pi/180)
#     end
# end

function disp_thickness_s(::TrippedN0012BoundaryLayer, alphastar)
    T = typeof(alphastar)
    # Equation 12 from the BPM report.
    alphastar_deg = alphastar*180/pi
    if alphastar_deg < 0
        # throw(DomainError(alphastar, "negative alphastar argument invalid"))
        return T(NaN)
    elseif alphastar_deg ≤ 5
        return 10^(0.0679*alphastar_deg)
    elseif alphastar_deg ≤ 12.5
        return 0.381*10^(0.1516*alphastar_deg)
    elseif alphastar_deg ≤ 25
        return 14.296*10^(0.0258*alphastar_deg)
    else
        # What should I do for angles of attack greater than 25°?
        # Maybe just keep the same thickness?
        return 14.296*10^(0.0258*25)*one(T)
    end
end

# function bl_thickness_s(::UntrippedN0012BoundaryLayer, alphastar)
#     # Equation 14 from the BPM report.
#     alphastar_deg = alphastar*180/pi
#     if alphastar_deg < 0
#         throw(DomainError(alphastar, "negative alphastar argument invalid"))
#     elseif alphastar_deg ≤ 7.5
#         return 10^(0.03114*alphastar)
#     elseif alphastar_deg ≤ 12.5
#         return 0.0303*10^(0.2336*alphastar)
#     elseif alphastar_deg ≤ 25
#         return 12*10^(0.0258*alphastar)
#     else
#         # What should I do for angles of attack greater than 25°?
#         # Maybe just keep the same thickness?
#         return 12*10^(0.0258*25*pi/180)
#     end
# end

function disp_thickness_s(::UntrippedN0012BoundaryLayer, alphastar)
    T = typeof(alphastar)
    # Equation 15 from the BPM report.
    alphastar_deg = alphastar*180/pi
    if alphastar_deg < 0
        # throw(DomainError(alphastar, "negative alphastar argument invalid"))
        return T(NaN)
    elseif alphastar_deg ≤ 7.5
        return 10^(0.0679*alphastar_deg)
    elseif alphastar_deg ≤ 12.5
        return 0.0162*10^(0.3066*alphastar_deg)
    elseif alphastar_deg ≤ 25
        return 52.42*10^(0.0258*alphastar_deg)
    else
        # What should I do for angles of attack greater than 25°?
        # Maybe just keep the same thickness?
        return 52.42*10^(0.0258*25)*one(T)
    end
end

function disp_thickness_top(bl::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}, alphastar)
    return ifelse(alphastar > alpha_zerolift(bl), disp_thickness_s(bl, alphastar), disp_thickness_p(bl, -alphastar))
end

function disp_thickness_bot(bl::Union{TrippedN0012BoundaryLayer,UntrippedN0012BoundaryLayer}, alphastar)
    return ifelse(alphastar > alpha_zerolift(bl), disp_thickness_p(bl, alphastar), disp_thickness_s(bl, -alphastar))
end

function bl_thickness_p(bl::AbstractBoundaryLayer, Re_c, alphastar)
    return bl_thickness_p(bl, alphastar)*bl_thickness_0(bl, Re_c)
end

# function bl_thickness_s(bl::AbstractBoundaryLayer, Re_c, alphastar)
#     return bl_thickness_s(bl, alphastar)*bl_thickness_0(bl, Re_c)
# end

function disp_thickness_bot(bl::AbstractBoundaryLayer, Re_c, alphastar)
    # (δ^*_p/δ^*_0)*(δ^*_0/c)
    return disp_thickness_bot(bl, alphastar)*disp_thickness_0(bl, Re_c)
end

function disp_thickness_top(bl::AbstractBoundaryLayer, Re_c, alphastar)
    return disp_thickness_top(bl, alphastar)*disp_thickness_0(bl, Re_c)
end
