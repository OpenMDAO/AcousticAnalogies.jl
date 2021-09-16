abstract type SourceElement end

@concrete struct CompactSourceElement <: SourceElement
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

    # Load *on the fluid* per unit span, and its time derivative.
    f0dot
    f1dot

    # Source time.
    τ

    # orientation of the element. Only used for WriteVTK.
    u
end

"""
    CompactSourceElement(ρ0, c0, r, θ, Δr, Λ, fn, fc, τ)

Construct a source element to be used with the compact form of Farassat's formulation 1A.

# Arguments
- ρ0: Ambient air density (kg/m^3)
- c0: Ambient speed of sound (m/s)
- r: radial coordinate of the element in the blade-fixed coordinate system (m)
- θ: angular offest of the element in the blade-fixed coordinate system (rad)
- Δr: length of the element (m)
- Λ: cross-sectional area of the element (m^2)
- fn: normal load per unit span *on the fluid* (N/m)
- fr: radial load *on the fluid* (N/m)
- fc: circumferential load *on the fluid* (N/m)
- τ: source time (s)
"""
function CompactSourceElement(ρ0, c0, r, θ, Δr, Λ, fn, fr, fc, τ)
    s, c = sincos(θ)
    y0dot = @SVector [0, r*c, r*s]
    T = eltype(y0dot)
    y1dot = @SVector zeros(T, 3)
    y2dot = @SVector zeros(T, 3)
    y3dot = @SVector zeros(T, 3)
    f0dot = @SVector [fn, c*fr - s*fc, s*fr + c*fc]
    T = eltype(f0dot)
    f1dot = @SVector zeros(T, 3)
    u = @SVector [0, c, s]

    return CompactSourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ, u)
end

"""
    (trans::KinematicTransformation)(se::CompactSourceElement)

Transform the motion and loading of a compact source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::CompactSourceElement)
    linear_only = false
    y0dot, y1dot, y2dot, y3dot = trans(se.τ, se.y0dot, se.y1dot, se.y2dot, se.y3dot, linear_only)
    linear_only = true
    f0dot, f1dot= trans(se.τ, se.f0dot, se.f1dot, linear_only)
    u = trans(se.τ, se.u, linear_only)

    return CompactSourceElement(se.ρ0, se.c0, se.Δr, se.Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, se.τ, u)
end

@concrete struct NonCompactSourceElement <: SourceElement
    # Ambient density.
    ρ0
    # Ambient speed of sound.
    c0
    # Area of surface element
    ΔA

    # Local density and its time derivative.
    ρ0dot
    ρ1dot

    # Normal unit vector and its time derivative.
    n0dot
    n1dot

    # Source position and its time derivatives.
    y0dot
    y1dot
    y2dot

    # Fluid velocity and its time derivative.
    u0dot
    u1dot

    # Load *on the fluid* per unit area, and its time derivative.
    f0dot
    f1dot

    # Source time.
    τ
end

"""
    (trans::KinematicTransformation)(se::NonCompactSourceElement)

Transform the motion and loading of a non-compact source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::NonCompactSourceElement)
    linear_only = false
    y0dot, y1dot, y2dot = trans(se.τ, se.y0dot, se.y1dot, se.y2dot, linear_only)
    y0dot, u0dot, u1dot = trans(se.τ, se.y0dot, se.u0dot, se.u1dot, linear_only)
    linear_only = true
    n0dot, n1dot = trans(se.τ, se.n0dot, n1dot)
    f0dot, f1dot= trans(se.τ, se.f0dot, se.f1dot, linear_only)

    return NonCompactSourceElement(se.ρ0, se.c0, se.ΔA, se.ρ0dot, se.ρ1dot, n0dot, n1dot, y0dot, y1dot, y2dot, u0dot, u1dot, f0dot, f1dot, se.τ)
end

@concrete struct StationaryNonCompactSourceElement <: SourceElement
    # Ambient density.
    ρ0
    # Ambient speed of sound.
    c0
    # Area of surface element
    ΔA

    # Local density and its time derivative.
    ρ0dot
    ρ1dot

    # Normal unit vector.
    n0dot

    # Source position.
    y0dot

    # Fluid velocity and its time derivative.
    u0dot
    u1dot

    # Load *on the fluid* per unit area, and its time derivative.
    f0dot
    f1dot

    # Source time.
    τ
end

"""
Supertype for an object that recieves a noise prediction when combined with an
acoustic analogy source; computational equivalent of a microphone.

    (obs::AcousticObserver)(t)

Calculate the position of the acoustic observer at time `t`.
"""
abstract type AcousticObserver end

"""
    StationaryAcousticObserver(x)

Construct an acoustic observer that does not move with position `x` (m).
"""
@concrete struct StationaryAcousticObserver <: AcousticObserver
    x
end

"""
    ConstVelocityAcousticObserver(t0, x0, v)

Construct an acoustic observer moving with a constant velocity `v`, located at
`x0` at time `t0`.
"""
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

"""
    adv_time(se::SourceElement, obs::AcousticObserver)

Calculate the time an acoustic wave emmitted by source `se` at time `se.τ` is
received by observer `obs`.
"""
adv_time(se::SourceElement, obs::AcousticObserver)

function adv_time(se::SourceElement, obs::StationaryAcousticObserver)
    rv = obs(se.τ) .- se.y0dot
    r = norm_cs_safe(rv)
    t = se.τ + r/se.c0
    return t
end

function adv_time(se::SourceElement, obs::ConstVelocityAcousticObserver)
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

"""
Output of the F1A calculation: the acoustic pressure value at time `t`, broken into monopole component `p_m` and
dipole component `p_d`.
"""
@concrete struct F1AOutput
    t
    p_m
    p_d
end


"""
    f1a(se::CompactSourceElement, obs::AcousticObserver, t_obs)

Calculate the acoustic pressure emitted by compact source element `se` and received by
observer `obs` at time `t_obs`, returning an [`F1AOutput`](@ref) object.

The correct value for `t_obs` can be found using [`adv_time`](@ref).
Alternatively you can omit the `t_obs` argument and it will be calculated.
"""
function f1a(se::CompactSourceElement, obs::AcousticObserver, t_obs)
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

    rhat1dot = -1/(r*r)*r1dot*rv + 1/r*rv1dot
    Mr1dot = (dot_cs_safe(rv2dot, rhat) + dot_cs_safe(rv1dot, rhat1dot))/(-se.c0)

    rhat2dot = (2/(r^3)*r1dot*r1dot*rv .- 1/(r^2)*r2dot*rv .- 2/(r^2)*r1dot*rv1dot .+ 1/r*rv2dot)

    Mr2dot = (dot_cs_safe(rv3dot, rhat) .+ 2*dot_cs_safe(rv2dot, rhat1dot) .+ dot_cs_safe(rv1dot, rhat2dot))/(-se.c0)

    # Rnm = r^(-n)*(1 - Mr)^(-m)
    R10 = 1/r
    R01 = 1/(1 - Mr)
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
    p_m = se.ρ0/(4*pi)*se.Λ*C1A*se.Δr

    # Dipole coefficients.
    D1A = R01*R11*rhat
    E1A = R01*(R11dot*rhat + R11*rhat1dot) + se.c0*R21*rhat

    # Dipole acoustic pressure!
    p_d = (dot_cs_safe(se.f1dot, D1A) + dot_cs_safe(se.f0dot, E1A))*se.Δr/(4*pi*se.c0)

    return F1AOutput(t_obs, p_m, p_d)
end

"""
    f1a(se::NonCompactSourceElement, obs::AcousticObserver, t_obs)

Calculate the acoustic pressure emitted by non-compact source element `se` and received by
observer `obs` at time `t_obs`, returning an [`F1AOutput`](@ref) object.

The correct value for `t_obs` can be found using [`adv_time`](@ref).
Alternatively you can omit the `t_obs` argument and it will be calculated.
"""
function f1a(se::NonCompactSourceElement, obs::AcousticObserver, t_obs)
    x_obs = obs(t_obs)

    rv = x_obs .- se.y0dot
    r = norm_cs_safe(rv)
    rhat = rv/r

    rv1dot = -se.y1dot
    r1dot = dot_cs_safe(rhat, rv1dot)

    rv2dot = -se.y2dot

    Mr = dot_cs_safe(-rv1dot/se.c0, rhat)

    rhat1dot = -1/(r*r)*r1dot*rv + 1/r*rv1dot
    Mr1dot = (dot_cs_safe(rv2dot, rhat) + dot_cs_safe(rv1dot, rhat1dot))/(-se.c0)

    # Rnm = r^(-n)*(1 - Mr)^(-m)
    R10 = 1/r
    R01 = 1/(1 - Mr)
    R11 = R10*R01
    R21 = R11*R10

    # Rnm1dot = d/dt(Rnm) = (-n*R10*r1dot + m*R01*Mr1dot)*Rnm
    R11dot = (-R10*r1dot + R01*Mr1dot)*R11

    A1 = R11
    A1dot = R11dot
    B1 = R11*rhat
    B1dot = R11dot*rhat + R11*rhat1dot
    C1 = se.c0*R21*rhat

    A1A = R01*A1
    B1A = R01*A1dot
    C1A = R01*B1
    D1A = R01*B1dot
    E1A = C1

    v = se.y1dot
    vdot = se.y2dot
    u = se.u0dot
    udot = se.u1dot
    n = se.n0dot
    ndot = se.n1dot
    v_n = dot_cs_safe(v, n)
    v_ndot = dot_cs_safe(vdot, n) + dot_cs_safe(v, ndot)
    u_n = dot_cs_safe(u, n)
    u_ndot = dot_cs_safe(udot, n) + dot_cs_safe(u, ndot)
    Q = se.ρ0*v_n + se.ρ0dot*(u_n - v_n)
    Qdot = se.ρ0*v_ndot + se.ρ1dot*(u_n - v_n) + se.ρ0dot*(u_ndot - v_ndot)
    p_m = (Qdot*A1A + Q*B1A)*(se.ΔA)/(4*pi)

    L = se.f0dot .+ se.ρ0dot*u*(u_n - v_n)
    Ldot = se.f1dot .+ se.ρ1dot*u*(u_n - v_n) .+ se.ρ0dot*udot*(u_n - v_n) .+ se.ρ0dot*u*(u_ndot - v_ndot)
    p_d = (dot_cs_safe(Ldot, C1A) + dot_cs_safe(L, D1A) + dot_cs_safe(L, E1A))*(se.ΔA)/(4*pi*se.c0)

    return F1AOutput(t_obs, p_m, p_d)
end

function f1a(se::StationaryNonCompactSourceElement, obs::AcousticObserver, t_obs)
    x_obs = obs(t_obs)

    rv = x_obs .- se.y0dot
    r = norm_cs_safe(rv)
    rhat = rv/r

    # Rnm = r^(-n)*(1 - Mr)^(-m)
    R10 = 1/r
    R01 = 1
    R11 = R10*R01
    R21 = R11*R10

    A1 = R11
    B1 = R11*rhat
    C1 = se.c0*R21*rhat

    A1A = R01*A1
    C1A = R01*B1
    E1A = C1

    u = se.u0dot
    udot = se.u1dot
    n = se.n0dot
    u_n = dot_cs_safe(u, n)
    u_ndot = dot_cs_safe(udot, n)
    Q = se.ρ0dot*u_n
    Qdot = se.ρ1dot*u_n + se.ρ0dot*u_ndot
    p_m = Qdot*A1A*(se.ΔA)/(4*pi)

    L = se.f0dot .+ se.ρ0dot*u*u_n
    Ldot = se.f1dot .+ se.ρ1dot*u*u_n .+ se.ρ0dot*udot*u_n .+ se.ρ0dot*u*u_ndot
    p_d = (dot_cs_safe(Ldot, C1A) + dot_cs_safe(L, E1A))*(se.ΔA)/(4*pi*se.c0)

    return F1AOutput(t_obs, p_m, p_d)
end

"""
    f1a(se::SourceElement, obs::AcousticObserver)

Calculate the acoustic pressure emitted by source element `se` and received by
observer `obs`, returning an [`F1AOutput`](@ref) object.
"""
function f1a(se::SourceElement, obs::AcousticObserver)
    t_obs = adv_time(se, obs)
    return f1a(se, obs, t_obs)
end

"""
    common_obs_time(apth::AbstractArray{<:F1AOutput}, period, n, axis=1)

Return a suitable time range for the collection of F1A acoustic pressures in `apth`.

The time range will begin near the latest start time of the acoustic pressures
in `apth`, and be an `AbstractVector` (really a `StepRangeLen`) of size `n` and
of time length `period`. `axis` indicates which axis of `apth` the time for a
source varies.
"""
function common_obs_time(apth, period, n, axis=1)
    # Make a single field struct array that behaves like a time array. 4%-6%
    # faster than creating the array with getproperty.
    t_obs = SingleFieldStructArray(apth, :t)

    # Get the first time for all the sources (returns a view ♥).
    t_starts = selectdim(t_obs, axis, 1)

    # Find the latest first time.
    t_common_start = ksmax(t_starts, 30/period)

    # Get the common observer time.
    dt = period/n
    t_common = t_common_start .+ (0:n-1)*dt

    return t_common
end

@concrete struct F1AAcousticPressure <: AcousticMetrics.AbstractAcousticPressure
    p_m
    p_d
    dt
    t0
end

"""
    F1AAcousticPressure([T=Float64,] n, dt, t0)

Construct an `F1AAcousticPressure` `struct` suitable for containing an acoustic prediction of length `n`, starting at time `t0` with time step `dt`.
"""
function F1AAcousticPressure(::Type{T}, n, dt, t0) where {T}
    p_m = Vector{T}(undef, n)
    p_d = Vector{T}(undef, n)
    return F1AAcousticPressure(p_m, p_d, dt, t0)
end

function F1AAcousticPressure(n, dt, t0)
    p_m = Vector{Float64}(undef, n)
    p_d = Vector{Float64}(undef, n)
    return F1AAcousticPressure(p_m, p_d, dt, t0)
end

"""
    F1AAcousticPressure(apth::AbstractArray{<:F1AOutput}, period::AbstractFloat, n::Integer, axis::Integer=1)

Construct an `F1AAcousticPressure` `struct` suitable for containing an acoustic prediction from an array of `F1AOutput` `struct`.

The elapsed time and length of the returned `F1AAcousticPressure` will be
`period` and `n`, respectively. `axis` indicates which axis the `apth` `struct`s
time varies. (`period`, `n`, `axis` are passed to [`common_obs_time`](@ref).)
"""
function F1AAcousticPressure(apth::AbstractArray{<:F1AOutput}, period, n, axis=1)
    # Get the common observer time.
    t_common = common_obs_time(apth, period, n, axis)

    # Allocate output arrays.
    T = typeof(first(apth).p_m)
    p_m = Vector{T}(undef, n)
    T = typeof(first(apth).p_d)
    p_d = Vector{T}(undef, n)

    # Create the output apth.
    dt = step(t_common)
    t0 = first(t_common)
    apth_out = F1AAcousticPressure(p_m, p_d, dt, t0)

    return apth_out
end

@inline AcousticMetrics.pressure(ap::F1AAcousticPressure) = ap.p_m + ap.p_d
@inline pressure_monopole(ap::F1AAcousticPressure) = ap.p_m
@inline pressure_dipole(ap::F1AAcousticPressure) = ap.p_d

"""
    combine!(apth_out::F1AAcousticPressure, apth::AbstractArray{<:F1AOutput}, axis; f_interp=akima)

Combine the acoustic pressures of multiple sources (`apth`) into a single acoustic pressure time history `apth_out`.

The input acoustic pressures `apth` are interpolated onto the time grid
returned by `time(apth_out)`. The interpolation is performed by the function `f_intep(xpt, ypt,
x)`, where `xpt` and `ytp` are the input grid and function values, respectively,
and `x` is the output grid.
"""
function combine!(apth_out, apth, axis; f_interp=akima)
    # This makes no difference compared to passing in a cache (an object with
    # working arrays that I'd copy stuff to) to this function (sometimes a
    # speedup of <1%, sometimes a slowdown of <1%). I'm sure it'd be worse if I
    # didn't pass in the cache. But it's nice to not have to worry about passing
    # it in.
    t_obs = SingleFieldStructArray(apth, :t)
    p_m = SingleFieldStructArray(apth, :p_m)
    p_d = SingleFieldStructArray(apth, :p_d)

    # Unpack the output arrays for clarity.
    t_common = AcousticMetrics.time(apth_out)
    p_m_interp = pressure_monopole(apth_out)
    p_d_interp = pressure_dipole(apth_out)

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

    return apth_out
end

"""
    combine(apth::AbstractArray{<:F1AOutput}, period::AbstractFloat, n::Integer, axis=1; f_interp=akima)

Combine the acoustic pressures of multiple sources (`apth`) into a single
acoustic pressure time history on a time grid of size `n` extending over time
length `period`.
"""
function combine(apth, period, n::Integer, axis::Integer=1; f_interp=akima)
    apth_out = F1AAcousticPressure(apth, period, n, axis)
    return combine!(apth_out, apth, axis; f_interp=f_interp)
end
