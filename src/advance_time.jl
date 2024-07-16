"""
    adv_time(se::AbstractCompactSourceElement, obs::AbstractAcousticObserver)

Calculate the time an acoustic wave emmited by source `se` at time `se.τ` is
recieved by observer `obs`.
"""
adv_time(se::AbstractCompactSourceElement, obs::AbstractAcousticObserver)

function adv_time(se::AbstractCompactSourceElement, obs::StationaryAcousticObserver)
    τ = source_time(se)
    rv = obs(τ) .- position(se)
    r = norm_cs_safe(rv)
    t = τ + r/speed_of_sound(se)
    return t
end

function adv_time(se::AbstractCompactSourceElement, obs::ConstVelocityAcousticObserver)
    # Source time of the source element.
    τ = source_time(se)

    # Ambient speed of sound of the source element.
    c0 = speed_of_sound(se)

    # Location of the observer at the source time.
    x = obs(τ)

    # Vector from the source to the observer at the source time.
    rv = x .- position(se)

    # Distance from the source to the observer at the source time.
    r = norm_cs_safe(rv)

    # Speed of the observer divided by speed of sound.
    Mo = norm_cs_safe(obs.v)/c0

    # Unit vector pointing from the source to the observer.
    rhat = rv/r

    # Velocity of observer dotted with rhat at the source time.
    Mor = dot_cs_safe(obs.v, rhat)/c0

    # Now get the observer time.
    t = τ + r/c0*((Mor + sqrt(Mor^2 + 1 - Mo^2))/(1 - Mo^2))

    return t
end

