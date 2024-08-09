abstract type AbstractDirectivity end
struct BPMDirectivity <: AbstractDirectivity end
struct BrooksBurleyDirectivity <: AbstractDirectivity end

abstract type AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler} <: AbstractCompactSourceElement end

"""
    doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver, t_obs)

Calculate the Doppler shift factor for noise emitted by source element `se` and recieved by observer `obs` at time `t_obs`, i.e. the ratio between an observer frequency `f` and emitted frequency `f_0`.

The correct value for `t_obs` can be found using [`adv_time`](@ref).
"""
function doppler_factor(se::AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,true}, obs::AbstractAcousticObserver, t_obs) where {TDirect,TUInduction,TMachCorrection}
    # Location of the observer at the observer time.
    x_obs = obs(t_obs)

    # Also need the speed of sound.
    c = speed_of_sound(se)

    # Get a unit vector pointing from the source position at the source time to the observer position at the observer time.
    rv = x_obs .- position(se)
    r = norm_cs_safe(rv)
    rhat = rv/r

    # So, now, if I dot the source velocity with `rhat`, that would give me the component of velocity of the source in the direction of the observer, positive if moving toward it, negative if moving away.
    v_src = dot_cs_safe(velocity(se), rhat)

    # And, if I dot the observer velocity `rhat`, that will give me the component of velocity of the observer in the direction of the source, positive if moving *away* from it, negative if moving toward.
    v_obs = dot_cs_safe(velocity(t_obs, obs), rhat)

    # Now we can get the factor.
    factor = (1 - v_obs/c) / (1 - v_src/c)

    return factor
end

function doppler_factor(se::AbstractBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,false}, obs::AbstractAcousticObserver, t_obs) where {TDirect,TUInduction,TMachCorrection}
    return 1
end

"""
    doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver)

Calculate the Doppler shift factor for noise emitted by source element `se` and recieved by observer `obs`, i.e. the ratio between an observer frequency `f` and emitted frequency `f_0`.

The correct value for `t_obs` will be found using [`adv_time`](@ref) internally.
"""
function doppler_factor(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver)
    # Do the advanced time calculation.
    t_obs = adv_time(se, obs)

    return doppler_factor(se, obs, t_obs)
end

function directivity(se::AbstractBroadbandSourceElement{BrooksBurleyDirectivity}, x_obs, top_is_suction)
    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    # But let's do it right, anyway!
    # if se.chord_cross_span_to_get_top_uvec
    #     # But, if the angle of attack is negative, then the "top" of the airfoil (which is normally the suction side) is actually the suction side.
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)
    #     else
    #         z_uvec_tmp = cross(se.span_uvec, se.chord_uvec)
    #     end
    # else
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.span_uvec, se.chord_uvec)
    #     else
    #         z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)
    #     end
    # end
    z_uvec_tmp = cross(se.chord_uvec, se.span_uvec)*ifelse(se.chord_cross_span_to_get_top_uvec, 1, -1)*ifelse(top_is_suction, 1, -1)
    z_uvec = z_uvec_tmp / norm_cs_safe(z_uvec_tmp)
    
    # Component of rv along the chord line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    x_er = dot_cs_safe(rv, se.chord_uvec)

    # Component of rv along the span line (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    y_er = dot_cs_safe(rv, se.span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find sin(Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(x_er/r_er)) = sqrt(r_er^2 - x_er^2)/r_er, and so sin(Θ_er)^2 = (r_er^2 - x_er^2)/r_er^2
    sin2Θer = (r_er^2 - x_er^2)/r_er^2

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 + z_er^2)

    # Need to find 2*sin(0.5*Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But there is a half-angle identity that says sin(θ/2)^2 = 0.5*(1 - cos(θ)).
    # So I actually want 2*sin(0.5*Θ_er)^2 = 2*0.5*(1 - cos(Θ_er)) = (1 - cos(Θ_er)).
    # But I can substitute in Θ_er = acos(x_er/r_er) and get 2*sin(0.5*Θ_er)^2 = 1 - x_er/r_er.
    twosin2halfΘer = 1 - x_er/r_er

    # Now just need the denominator: (1 - M_tot*cos(ξR))^4.
    # M_tot is the "total" velocity from... hmm... what perspective?
    # Let's see... it looks like it's suppose to be from the fluid, aka the global frame.
    # The definition is Brooks and Burley AIAA 2001-2210, equation (14):
    #
    #   V_tot = V - V_wt - V_ind
    #
    # where 
    #
    #   * V is the velocity due to the rotation of the blade element
    #   * V_wt is the wind tunnel velocity, which is positive when it goes against the motion of the blade element.
    #   * V_ind is "the induced velocity due to the near and far wake of the rotor," and appears to be positive in roughly the thrust direction.
    #
    # So if I calculate V - V_wt, that's the "actual" velocity of the blade element, i.e., the velocity of the blade element relative to the fluid far away from the blade element, since it doesn't include the induced velocity.
    # That's what I usually think of as the "actual" velocity, since it's what a stationary observer would observe on a calm day.
    # But when we add in the induced velocity, I think what we're finding is the velocity of the blade element relative to the nearfield velocity.
    # Cool.
    # But does that mean I add or subtract `se.y1dot_fluid` from `se.y1dot`?
    # Well, let's think about that.
    # First, let's say I start with stuff in the blade-fixed frame.
    # And let's say I'm imagining that, from the global frame, I'm assuming the
    # blade element is moving in the positive x direction, initially aligned
    # with the y axis
    # So I think all I need to do is just use se.y1dot_fluid + se.y1dot.
    # Now, cos(ξ_r) is defined by equation (18) in Brooks and Burley AIAA 2001-2210, which is the angle between the radiation vector (rv here) and the total velocity (se.y1dot here).
    # But I can simplify that by just finding the unit radiation vector, then dotting that with the velocity vector, and dividing by c0.

    # Unit radiation vector.
    r_uvec = rv./r_er
    # Equation 14 from Brooks and Burley AIAA 2001-2210.
    Vtotal = se.y1dot - se.y1dot_fluid
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Convective amplification factor for the two directivity functions.
    conv_amp = 1/(1 - Mtotcosξr)^4

    # Now I can finally find the directivity function!
    # Equation (19) from Brooks and Burley AIAA 2001-2210.
    # Dl = (sin2Θer*sin2Φer)/(1 - Mtotcosξr)^4
    Dl = (sin2Θer*sin2Φer)*conv_amp

    # Now I can finally find the directivity function!
    # Equation (20) from Brooks and Burley AIAA 2001-2210.
    # Dh = (twosin2halfΘer*sin2Φer)/(1 - Mtotcosξr)^4
    Dh = (twosin2halfΘer*sin2Φer)*conv_amp

    return r_er, Dl, Dh
end

function directivity(se::AbstractBroadbandSourceElement{BPMDirectivity}, x_obs, top_is_suction)
    # Position vector from source to observer.
    rv = x_obs .- se.y0dot

    # Distance from source to observer.
    r_er = norm_cs_safe(rv)

    # So, the BPM report uses the local flow velocity, not the chord line, to define the x direction.
    # So, I want to get a unit vector in that direction.
    # Should it include induction?
    # It won't matter for comparing to the data in the BPM report, since the flow including and excluding induction would be in the same direction.
    # In the BPM report Appendix B, the x direction is defined as the opposite of the motion of the source element/flat plate, so I guess I won't use induction.
    # But I want the velocity to be normal to the span direction, so let's remove_that.
    # So, want the x direction to be opposite the velocity of the source element.
    x_vec_tmp1 = -se.y1dot
    # Then we want to remove any part of the velocity in the direction of the span.
    x_vec_tmp2 = x_vec_tmp1 - dot_cs_safe(x_vec_tmp1, se.span_uvec)*x_vec_tmp1
    # Now make it a unit vector:
    x_uvec = x_vec_tmp2 / norm_cs_safe(x_vec_tmp2)

    # Unit vector normal to both the span and chord directions.
    # Does the order matter?
    # Doesn't look like it, since we're only using it to find z_er, which we square.
    # But it's supposed to be pointing from the pressure to the suction side, which we can figure out, so let's do it the right way.
    # if se.chord_cross_span_to_get_top_uvec
    #     if top_is_suction
    #         z_uvec_tmp = cross(x_uvec, se.span_uvec)
    #     else
    #         z_uvec_tmp = cross(se.span_uvec, x_uvec)
    #     end
    # else
    #     if top_is_suction
    #         z_uvec_tmp = cross(se.span_uvec, x_uvec)
    #     else
    #         z_uvec_tmp = cross(x_uvec, se.span_uvec)
    #     end
    # end
    z_uvec_tmp = cross(x_uvec, se.span_uvec)*ifelse(se.chord_cross_span_to_get_top_uvec, 1, -1)*ifelse(top_is_suction, 1, -1)
    z_uvec = z_uvec_tmp / norm_cs_safe(z_uvec_tmp)
    
    # Component of rv along the chord line (see Figure B3 in the BPM report).
    x_er = dot_cs_safe(rv, x_uvec)

    # Component of rv along the span line (see Figure B3 the BPM report).
    y_er = dot_cs_safe(rv, se.span_uvec)

    # Component of rv in the direction normal to both span and chord (see Figure 11 in Brooks and Burley AIAA 2001-2210).
    z_er = dot_cs_safe(rv, z_uvec)

    # Need to find sin(Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(x_er/r_er)) = sqrt(r_er^2 - x_er^2)/r_er, and so sin(Θ_er)^2 = (r_er^2 - x_er^2)/r_er^2
    sin2Θer = (r_er^2 - x_er^2)/r_er^2

    # Need to find sin(Φ_er)^2, where Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2)), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But sin(acos(y_er/sqrt(y_er^2 + z_er^2))) = z_er/sqrt(y_er^2 + z_er^2), and so sin(Φ_er)^2 = z_er^2/(y_er_^2 + z_er^2).
    sin2Φer = (z_er^2)/(y_er^2 + z_er^2)

    # Need to find 2*sin(0.5*Θ_er)^2, where Θ_er = acos(x_er/r_er), equation (21) from Brooks and Burley AIAA 2001-2210.
    # But there is a half-angle identity that says sin(θ/2)^2 = 0.5*(1 - cos(θ)).
    # So I actually want 2*sin(0.5*Θ_er)^2 = 2*0.5*(1 - cos(Θ_er)) = (1 - cos(Θ_er)).
    # But I can substitute in Θ_er = acos(x_er/r_er) and get 2*sin(0.5*Θ_er)^2 = 1 - x_er/r_er.
    twosin2halfΘer = 1 - x_er/r_er

    # Now just need the denominator: (1 - M_tot*cos(ξR))^4.
    # M_tot is the "total" velocity from... hmm... what perspective?
    # Let's see... it looks like it's suppose to be from the fluid, aka the global frame.
    # The definition is Brooks and Burley AIAA 2001-2210, equation (14):
    #
    #   V_tot = V - V_wt - V_ind
    #
    # where 
    #
    #   * V is the velocity due to the rotation of the blade element
    #   * V_wt is the wind tunnel velocity, which is positive when it goes against the motion of the blade element.
    #   * V_ind is "the induced velocity due to the near and far wake of the rotor," and appears to be positive in roughly the thrust direction.
    #
    # So if I calculate V - V_wt, that's the "actual" velocity of the blade element, i.e., the velocity of the blade element relative to the fluid far away from the blade element, since it doesn't include the induced velocity.
    # That's what I usually think of as the "actual" velocity, since it's what a stationary observer would observe on a calm day.
    # But when we add in the induced velocity, I think what we're finding is the velocity of the blade element relative to the nearfield velocity.
    # Cool.
    # But does that mean I add or subtract `se.y1dot_fluid` from `se.y1dot`?
    # Well, let's think about that.
    # First, let's say I start with stuff in the blade-fixed frame.
    # And let's say I'm imagining that, from the global frame, I'm assuming the
    # blade element is moving in the positive x direction, initially aligned
    # with the y axis
    # So I think all I need to do is just use se.y1dot_fluid + se.y1dot.
    # Now, cos(ξ_r) is defined by equation (18) in Brooks and Burley AIAA 2001-2210, which is the angle between the radiation vector (rv here) and the total velocity (se.y1dot here).
    # But I can simplify that by just finding the unit radiation vector, then dotting that with the velocity vector, and dividing by c0.

    # Unit radiation vector.
    r_uvec = rv./r_er
    # For the BPM directivity function, the velocity/Mach number doesn't include induction.
    Vtotal = se.y1dot
    # Mach number vectory in the direction of the radiation vector.
    Mtotcosξr = dot_cs_safe(Vtotal, r_uvec)/se.c0

    # Convective amplification factor for the low-freqency directivity function.
    conv_amp_l = 1/(1 - Mtotcosξr)^4

    # The BPM high-frequency convective amplification factor is a bit different.
    # It has a factor (M - M_c)*cos(Θ_er), which, in the more general coordinate system of Brooks & Burley would be -(M - M_c)*cos(ξ_r).
    # So, the `M` is the speed of the blade element without induction, and `M_c` is the velocity of the blade element including induction.
    # So `M_c = se.y1dot - se.y1dot_fluid` and then `M - M_c = se.y1dot - (se.y1dot - se.y1dot_fluid) = se.y1dot_fluid.
    # And so what we'd want to do is this:
    M_minus_M_ccosξr = dot_cs_safe(se.y1dot_fluid, r_uvec)/se.c0
    conv_amp_h = 1/((1 - Mtotcosξr)*(1 - M_minus_M_ccosξr)^2)

    # Now I can finally find the directivity function!
    # Equation (B2) from the BPM report.
    Dl = (sin2Θer*sin2Φer)*conv_amp_l

    # Now I can finally find the directivity function!
    # Equation (B1) from the BPM report.
    Dh = (twosin2halfΘer*sin2Φer)*conv_amp_h

    return r_er, Dl, Dh
end

function angle_of_attack(se::AbstractBroadbandSourceElement)
    # Find the total velocity of the fluid from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    # Vtotal = -(se.y1dot - se.y1dot_fluid)
    Vtotal = se.y1dot_fluid - se.y1dot 

    # To get the angle of attack, I need to find the components of the velocity in the chordwise direction, and the direction normal to both the chord and span.
    # So, first need to get a vector normal to both the chord and span, pointing from pressure side to suction side.
    normal_uvec_tmp = ifelse(se.chord_cross_span_to_get_top_uvec,
        cross(se.chord_uvec, se.span_uvec),
        cross(se.span_uvec, se.chord_uvec))
    normal_uvec = normal_uvec_tmp ./ norm_cs_safe(normal_uvec_tmp)

    # Now get the component of velocity in the chord_uvec and normal_uvec directions.
    V_chordwise = dot_cs_safe(Vtotal, se.chord_uvec)
    V_normal = dot_cs_safe(Vtotal, normal_uvec)

    # Now we can find the angle of attack.
    alphastar = atan_cs_safe(V_normal, V_chordwise)
    # alphastar = atan(V_normal, V_chordwise)
    
    return alphastar
end

function speed_normal_to_span(se::AbstractBroadbandSourceElement{TDirect,true}) where {TDirect}
    # Find the total velocity of the fluid including induction, from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    Vtotal = se.y1dot_fluid - se.y1dot
    # Find the component of the velocity in the direction of the span.
    Vspan = dot_cs_safe(Vtotal, se.span_uvec)*se.span_uvec
    # Subtract that from the total velocity to get the velocity normal to the span, then get the norm for the speed normal to span.
    return norm_cs_safe(Vtotal - Vspan)
end

function speed_normal_to_span(se::AbstractBroadbandSourceElement{TDirect,false}) where {TDirect}
    # Find the total velocity of the fluid, not including induction, from the perspective of the blade element, which is just the total velocity of the blade element from the perspective of the fluid with the sign switched.
    Vtotal = -se.y1dot
    # Find the component of the velocity in the direction of the span.
    Vspan = dot_cs_safe(Vtotal, se.span_uvec)*se.span_uvec
    # Subtract that from the total velocity to get the velocity normal to the span, then get the norm for the speed normal to span.
    return norm_cs_safe(Vtotal - Vspan)
end

function noise(se::AbstractBroadbandSourceElement, obs::AbstractAcousticObserver, freqs::AcousticMetrics.AbstractProportionalBands{3, :center})
    t_obs = adv_time(se, obs)
    return noise(se, obs, t_obs, freqs)
end

