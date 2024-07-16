"""
Supertype for an object that recieves a noise prediction when combined with an
acoustic analogy source; computational equivalent of a microphone.

    (obs::AbstractAcousticObserver)(t)

Calculate the position of the acoustic observer at time `t`.
"""
abstract type AbstractAcousticObserver end

"""
    StationaryAcousticObserver(x)

Construct an acoustic observer that does not move with position `x` (m).
"""
struct StationaryAcousticObserver{Tx} <: AbstractAcousticObserver
    x::Tx
end

"""
    velocity(t_obs, obs::StationaryAcousticObserver)

Return the velocity of `obs` at time `t_obs` (hint—will always be zero ☺)
"""
@inline velocity(t_obs, obs::StationaryAcousticObserver) = zero(obs.x)

"""
    ConstVelocityAcousticObserver(t0, x0, v)

Construct an acoustic observer moving with a constant velocity `v`, located at
`x0` at time `t0`.
"""
struct ConstVelocityAcousticObserver{Tt0,Tx0,Tv} <: AbstractAcousticObserver
    t0 ::Tt0
    x0::Tx0
    v::Tv
end

function (obs::StationaryAcousticObserver)(t)
    return obs.x
end

function (obs::ConstVelocityAcousticObserver)(t)
    return obs.x0 .+ (t - obs.t0).*obs.v
end

"""
    velocity(t_obs, obs::ConstVelocityAcousticObserver)

Return the velocity of `obs` at time `t_obs` (hint—will always be the same)
"""
@inline velocity(t_obs, obs::ConstVelocityAcousticObserver) = obs.v
