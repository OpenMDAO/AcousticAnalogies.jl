module AcousticAnalogies

using ConcreteStructs: @concrete
using FLOWMath: akima, linear
using KinematicCoordinateTransformations
using LinearAlgebra: cross, norm, mul!
using StaticArrays

export CompactSourceElement
export AcousticObserver, StationaryAcousticObserver, ConstVelocityAcousticObserver
export AcousticPressure
export adv_time
export create_cache, f1a
export common_obs_time!, common_obs_time
export combine!, combine

include("utils.jl")

@concrete struct CompactSourceElement
    # Density.
    ρ0
    # Speed of sound.
    c0
    # Radial length of element.
    Δr
    # Cross-sectional area.
    Λ
    # Source position and its time derivatives.
    y0dot
    y1dot
    y2dot
    y3dot

    # Load *on the fluid*, and its time derivative.
    f0dot
    f1dot

    # Source time.
    τ
end

function CompactSourceElement(ρ0, c0, r, θ, Δr, Λ, fn, fc, τ)
    T = typeof(r)
    y0dot = @SVector [zero(T), r*cos(θ), r*sin(θ)]
    y1dot = @SVector zeros(T, 3)
    y2dot = @SVector zeros(T, 3)
    y3dot = @SVector zeros(T, 3)
    T = typeof(fn)
    f0dot = @SVector [-fn, -sin(θ)*fc, cos(θ)*fc]
    f1dot = @SVector zeros(T, 3)

    return CompactSourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ)
end

function CompactSourceElement(ρ0, c0, r, Δr, Λ, fn, fc, τ)
    θ = zero(r)
    return CompactSourceElement(ρ0, c0, r, θ, Δr, Λ, fn, fc, τ)
end

function (trans::KinematicTransformation)(se::CompactSourceElement)
    linear_only = false
    y0dot, y1dot, y2dot, y3dot = trans(se.τ, se.y0dot, se.y1dot, se.y2dot, se.y3dot, linear_only)
    linear_only = true
    f0dot, f1dot= trans(se.τ, se.f0dot, se.f1dot, linear_only)

    return CompactSourceElement(se.ρ0, se.c0, se.Δr, se.Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, se.τ)
end

abstract type AcousticObserver end

@concrete struct StationaryAcousticObserver <: AcousticObserver
    x
end

@concrete struct ConstVelocityAcousticObserver <: AcousticObserver
    t0 
    x0
    v
end

function (obs::StationaryAcousticObserver)(t)
    return obs.x
end

function (obs::ConstVelocityAcousticObserver)(t)
    return obs.x0 .+ (t - obs.t0).*obs.v
end

function adv_time(se::CompactSourceElement, obs::StationaryAcousticObserver)
    rv = obs.x .- se.y0dot
    r = norm_cs_safe(rv)
    t = se.τ + r/se.c0
    return t
end

function adv_time(se::CompactSourceElement, obs::ConstVelocityAcousticObserver)
    # Location of the observer at the source time.
    x = obs(se.τ)

    # Vector from the source to the observer at the source time.
    rv = x .- se.y0dot

    # Distance from the source to the observer at the source time.
    r = norm_cs_safe(rv)

    # Speed of the observer divided by speed of sound.
    Mo = norm_cs_safe(obs.v)/se.c0

    # Unit vector pointing from the source to the observer.
    rhat = rv/r

    # Velocity of observer dotted with rhat at the source time.
    Mor = dot_cs_safe(obs.v, rhat)/se.c0

    # Now get the observer time.
    t = se.τ + r/se.c0*((Mor + sqrt(Mor^2 + 1 - Mo^2))/(1 - Mo^2))

    return t
end

@concrete struct AcousticPressure
    t
    p_m
    p_d
end

function create_cache(apth::AbstractArray{AcousticPressure{T1,T2,T3}, N}) where {T1,T2,T3,N}
    sz = size(apth)
    t = Array{T1}(undef, sz)
    p_m = Array{T2}(undef, sz)
    p_d = Array{T3}(undef, sz)
    return AcousticPressure(t, p_m, p_d)
end

function f1a(se::CompactSourceElement, obs, t_obs)
    x_obs = obs(t_obs)

    rv = x_obs .- se.y0dot
    r = norm_cs_safe(rv)
    rhat = rv/r

    rv1dot = -se.y1dot
    r1dot = dot_cs_safe(rhat, rv1dot)

    rv2dot = -se.y2dot
    r2dot = (dot_cs_safe(rv1dot, rv1dot) + dot_cs_safe(rv, rv2dot) - r1dot*r1dot)/r

    rv3dot = -se.y3dot

    Mr = dot_cs_safe(-rv1dot/se.c0, rhat)

    rhat1dot = -1.0/(r*r)*r1dot*rv + 1.0/r*rv1dot
    Mr1dot = (dot_cs_safe(rv2dot, rhat) + dot_cs_safe(rv1dot, rhat1dot))/(-se.c0)

    rhat2dot = (2.0/(r^3)*r1dot*r1dot*rv .- 1.0/(r^2)*r2dot*rv .- 2.0/(r^2)*r1dot*rv1dot .+ 1.0/r*rv2dot)

    Mr2dot = (dot_cs_safe(rv3dot, rhat) .+ 2.0*dot_cs_safe(rv2dot, rhat1dot) .+ dot_cs_safe(rv1dot, rhat2dot))/(-se.c0)

    # Rnm = r^(-n)*(1.0 - Mr)^(-m)
    R10 = 1.0/r
    R01 = 1.0/(1.0 - Mr)
    R11 = R10*R01
    R02 = R01*R01
    R21 = R11*R10

    # Rnm1dot = d/dt(Rnm) = (-n*R10*r1dot + m*R01*Mr1dot)*Rnm
    R10dot = -R10*r1dot*R10
    R01dot = R01*Mr1dot*R01
    R11dot = (-R10*r1dot + R01*Mr1dot)*R11

    R11dotdot = (-R10dot*r1dot - R10*r2dot + R01dot*Mr1dot + R01*Mr2dot)*R11 + (-R10*r1dot + R01*Mr1dot)*R11dot

    # Monopole coefficient.
    C1A = R02*R11dotdot + R01*R01dot*R11dot

    # Monople acoustic pressure!
    p_m = se.ρ0/(4.0*pi)*se.Λ*C1A*se.Δr

    # Dipole coefficients.
    D1A = R01*R11*rhat
    E1A = R01*(R11dot*rhat + R11*rhat1dot) + se.c0*R21*rhat

    # Dipole acoustic pressure!
    p_d = (dot_cs_safe(se.f1dot, D1A) + dot_cs_safe(se.f0dot, E1A))*se.Δr/(4.0*pi*se.c0)

    return AcousticPressure(t_obs, p_m, p_d)
end

function f1a(se::CompactSourceElement, obs::T) where {T<:AcousticObserver}
    t_obs = adv_time(se, obs)
    return f1a(se, obs, t_obs)
end

function common_obs_time!(t_common, t_obs, period, axis=1)
    n = length(t_common)

    # Get the first time for all the sources (returns a view ♥).
    t_starts = selectdim(t_obs, axis, 1)

    # Get one of the start times for shifting (makes the smooth_max function's
    # numerical behavior better).
    t0 = t_starts[1]

    # Find the latest first time.
    t_common_start = smooth_max(t_starts .- t0, k=30.0/period) + t0

    # Get the common observer time.
    dt = period/n
    t_common .= t_common_start .+ (0:n-1)*dt

    return nothing
end

function common_obs_time!(t_common, apth::AbstractArray{<:AcousticPressure}, period, axis=1)
    n = length(t_common)

    # Get the first acoustic pressure for all the sources (returns a view ♥).
    apth_starts = selectdim(apth, axis, 1)

    # Get one of the start times for shifting (makes the smooth_max function's
    # numerical behavior better).
    t0 = apth_starts[1].t

    # Find the latest first time.
    t_common_start = smooth_max_time(apth_starts, 30.0/period, t0) + t0

    # Get the common observer time.
    dt = period/n
    t_common .= t_common_start .+ (0:n-1)*dt

    return nothing
end

function common_obs_time(apth, period, n, axis=1)
    T = typeof(first(apth).t)
    t_common = Vector{T}(undef, n)

    common_obs_time!(t_common, apth, period, axis)

    return t_common
end

function combine!(apth_out, apth, axis::Integer=1; f_interp=akima, cache=create_cache(apth))

    # Set the values for the cache.
    for I in eachindex(apth)
        cache.t[I] = apth[I].t
        cache.p_m[I] = apth[I].p_m
        cache.p_d[I] = apth[I].p_d
    end

    # Unpack the cache for clarity.
    t_obs = cache.t
    p_m = cache.p_m
    p_d = cache.p_d

    # Unpack the output arrays for clarity.
    t_common = apth_out.t
    p_m_interp = apth_out.p_m
    p_d_interp = apth_out.p_d

    dimsAPTH = [axes(t_obs)...]
    ndimsAPTH = ndims(t_obs)
    alldims = [1:ndimsAPTH;]  # Is this any better than `collect(1:ndimsAPTH)`?

    otherdims = setdiff(alldims, axis)
    itershape = tuple(dimsAPTH[otherdims]...)

    idx = Any[first(ind) for ind in axes(t_obs)]
    idx[axis] = Colon()

    nidx = length(otherdims)
    indices = CartesianIndices(itershape)

    # Zero out the output arrays.
    fill!(p_m_interp, zero(eltype(p_m_interp)))
    fill!(p_d_interp, zero(eltype(p_d_interp)))

    # Loop through the indices.
    for I in indices
        for i in 1:nidx
            idx[otherdims[i]] = I.I[i]
        end
        # Now I have the current indices of the source that I want to
        # interpolate.
        p_m_interp .+= f_interp(t_obs[idx...], p_m[idx...], t_common)
        p_d_interp .+= f_interp(t_obs[idx...], p_d[idx...], t_common)
    end

    return nothing
end

function combine(apth, t_common::AbstractArray, axis::Integer=1; f_interp=akima)
    # Allocate output arrays.
    nout = length(t_common)
    T = typeof(first(apth).p_m)
    p_m_interp = zeros(T, nout)
    T = typeof(first(apth).p_d)
    p_d_interp = zeros(T, nout)

    # Create the output apth.
    apth_out = AcousticPressure(t_common, p_m_interp, p_d_interp)

    # Do it.
    combine!(apth_out, apth, axis; f_interp=f_interp)

    return apth_out
end

function combine(apth, period::AbstractFloat, n::Integer, axis::Integer=1; f_interp=akima)
    # Get a common time grid.
    t_common = common_obs_time(apth, period, n, axis)
    return combine(apth, t_common, axis; f_interp=f_interp)
end
end # module
