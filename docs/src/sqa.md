```@meta
CurrentModule = AADocs
```
# Software Quality Assurance

## Tests
AcousticAnalogies.jl uses the usual Julia testing framework to implement and run tests.
The tests can be run locally after installing AcousticAnalogies.jl, and are also run automatically on GitHub Actions.

To run the tests locally, from the Julia REPL, type `]` to enter the Pkg prompt, then

```julia-repl
(jl_jncZ1E) pkg> test AcousticAnalogies
     Testing Running tests...
Test Summary:       | Pass  Total  Time
Advanced time tests |    2      2  7.0s
Test Summary:           | Pass  Total  Time
Combine F1AOutput tests |    8      8  2.7s
Test Summary: | Pass  Total  Time
F1A tests     |    2      2  5.9s
Test Summary:               | Pass  Total  Time
CCBlade private utils tests |    1      1  0.3s
Test Summary:                     | Pass  Total  Time
CCBlade CompactSourceElement test |   12     12  3.4s
Test Summary:     | Pass  Total  Time
ANOPP2 Comparison |  176    176  5.9s
Test Summary:    | Time
ForwardDiff test | None  14.1s
     Testing AcousticAnalogies tests passed

(jl_jncZ1E) pkg> 
```

(The output associated with installing all the dependencies the tests need aren't shown above.)

Here is a description of each category of test:

### Advanced Time Tests
The F1A calculation is concerned with roughly two types of objects: acoustic sources and acoustic observers.
Acoustic sources are things that make noise, and, for AcousticAnalogies.jl, would typically be a portion of some type of aerodynamic lifting surface (like a propeller blade).
An acoustic observer is just a fancy name for a person or microphone that will hear the noise emitted by the source.
Both the source and observer may be stationary, but more likely will be moving.

During the F1A calculation, we need to know the time at which an acoustic wave emitted by the source encounters the observer.
Mathematically, we need to solve the equation

```math
R(t) = t - \left( \tau + \frac{|\vec{x}(t) - \vec{y}(τ)|}{c_0} \right) = 0
```

whereon

  * ``τ`` is the time the source has emitted an acoustic disturbance
  * ``t`` is the time the observer encounters the acoustic disturbance
  * ``\vec{y}`` is the position of the source
  * ``\vec{x}`` is the position of the observer
  * ``c_0`` is the speed of sound

AcousticAnalogies.jl currently uses an advanced time approach to solving this equation.
This means we start with knowledge of ``\tau`` and then calculate ``t``—we "advance" the source time to the observer time by adding the amount of time it takes for the acoustic disturbance to travel from ``y`` to ``x``.

Now, the ``R(t) = 0`` equation is quite easy to solve if the observer is stationary.
In that case, ``x`` is not a function of ``t``, and so solving for ``t`` just involves moving everything in the parenthesis to the right-hand side.
But if the observer is moving, things are more complicated.
It may be impossible to solve for ``t`` explicitly in that case.
It turns out, however, that there is an explicit solution for ``t`` in the advanced time approach if the observer is moving at a constant rate (see D. Casolino [http://dx.doi.org/10.1016/S0022-460X(02)00986-0](http://dx.doi.org/10.1016/S0022-460X(02)00986-0)).
The constant velocity case is actually quite handy, since it's what we need to compare to wind tunnel data.

So, how do we test that we've implemented the solution to the ``R(t) = 0`` advanced time equation correctly?
In AcousticAnalogies.jl, we just use the nonlinear solver provided by [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl), and compare its solution to AcousticAnalogies.jl.
Here's how to do that:

```@example adv_time_tests
using AcousticAnalogies: AcousticAnalogies
using LinearAlgebra: norm
using NLsolve: NLsolve
using StaticArrays

# Create a source element for the test.
# The only things about the source element that matters to the advanced # time calculation is the time and position, and the speed of sound.
# So everything else will be take on dummy values.
τ = 2.5
y = @SVector [-4.0, 3.0, 6.0]
c0 = 2.0
dummy0 = 1.0
dummy3 = @SVector [0.0, 0.0, 0.0]
se = AcousticAnalogies.CompactSourceElement(dummy0, c0, dummy0, dummy0, y, dummy3, dummy3, dummy3, dummy3, dummy3, τ, dummy3)

# Define a function that solves the advanced time equation using `nlsolve.
function adv_time_nlsolve(se, obs)
    # Create the residual equation that we'll solve.
    # nlsolve assumes the residual function takes in and returns arrays.
    R(t) = [t[1] - (se.τ + norm(obs(t[1]) .- se.y0dot)/se.c0)]

    # Solve the advanced time equation.
    result = NLsolve.nlsolve(R, [1.0], autodiff=:forward)
    if !NLsolve.converged(result)
        @error "nlsolve advanced time calculation did not converge:\n$(result)"
    end
    t_obs = result.zero[1]
    return t_obs
end

# Let's try it out.

# First, a stationary observer:
x0 = @SVector [-3.0, 2.0, 8.5]
obs = AcousticAnalogies.StationaryAcousticObserver(x0)
t_exact = AcousticAnalogies.adv_time(se, obs)
t_nlsolve = adv_time_nlsolve(se, obs)
println("stationary observer, exact: $(t_exact), nlsorve: $(t_nlsolve), difference = $(t_exact - t_nlsolve)")

# Next, a constant velocity observer:
t0 = 3.5
x0 = @SVector [-2.0, 3.5, 6.25]
v = @SVector [-1.5, 1.5, 3.25]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0, x0, v)
t_exact = AcousticAnalogies.adv_time(se, obs)
t_nlsolve = adv_time_nlsolve(se, obs)
println("constant velocity observer, exact: $(t_exact), nlsorve: $(t_nlsolve), difference = $(t_exact - t_nlsolve)")
```

Almost identical results, so things are good!

### Combine `F1AOutput` Tests
The function `f1a(se::CompactSourceElement, obs::AcousticObserver)` uses Farassat's formulation 1A to perform a prediction of the noise experienced by one observer `obs` due to one acoustic source `se`.
Typically we will not have just one source, however.
For example, the [guided example in the docs](@ref guided_example) uses 30 "source elements" to model each propeller blade.
But we're interested in the acoustics experienced by `obs` due to **all** of the source elements, not just one.
So, we need to combine the output of `f1a` for one observer and all of the source elements.
In AcousticAnalogies.jl, this is done by interpolating the time history of each source element's acoustics (the "pressure time history") onto a common chunk of time, and then adding them up.
No big deal.

But, how do we test the "interpolating and adding" routine, aka [`AcousticAnalogies.combine`](@ref)?
That's pretty simple, actually: we just define some arbitrary functions that we'll use to create some pressure time histories, add them using the `combine` routine, and then compare that result to those created via evaluating those arbitrary functions on the same time grid used by the `combine` routine.
If those match, then the test passes, and everything in the `combine` routine should be good.
Let's try that:

```@example combine_test
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using GLMakie
using Random

# Goal is to verify that the code can faithfully combine two acoustic pressures on different time "grids" onto a single common grid.
# These will be our made up functions:
fa(t) = sin(2*pi*t) + 0.2*cos(4*pi*(t-0.1))
fb(t) = cos(6*pi*t) + 0.3*sin(8*pi*(t-0.2))

# Now we'll make some made up time grids.
n = 101
t1 = collect(range(0.0, 1.0, length=n))
dt = t1[2] - t1[1]
# Add a bit of random noise to the time grid.
# Make sure that the amount of # randomness isn't large enough to make the time values non-monotonically increasing (i.e., they don't overlap).
noise = 0.49.*dt.*(1 .- 2 .* rand(size(t1)...))
t1 .+= noise

t2 = collect(range(0.1, 1.1, length=n))
dt = t2[2] - t2[1]
t2 .+= 0.49.*dt.*(1 .- 2 .* rand(size(t2)...))

# Now let's create a bunch of pressure time histories on the time grids we just defined.
apth1 = @. AcousticAnalogies.F1AOutput(t1, fa(t1), 2*fa(t1))
apth2 = @. AcousticAnalogies.F1AOutput(t2, fb(t2), 3*fb(t2))

# Calculate the "exact" answer by coming up with a common time, then evaluating the test functions directly on the common time grid.
period = 0.5
n_out = 51
t_start = max(t1[1], t2[1])
t_common = t_start .+ (0:n_out-1).*(period/n_out)

p_m = @. fa(t_common)+fb(t_common)
p_d = @. 2*fa(t_common)+3*fb(t_common)

even_length = iseven(n_out)
apth_test = AcousticAnalogies.F1APressureTimeHistory{even_length}(p_m, p_d, step(t_common), first(t_common))

# Put all the acoustic pressures in one array.
apth = hcat(apth1, apth2)

# Combine.
apth_out = AcousticAnalogies.combine(apth, period, n_out)

# Plot the two solutions.
fig2 = Figure()
ax2_1 = fig2[1, 1] = Axis(fig2, xlabel="time", ylabel="acoustic pressure, monopole")
ax2_2 = fig2[2, 1] = Axis(fig2, xlabel="time", ylabel="acoustic pressure, dipole")
scatter!(ax2_1, AcousticMetrics.time(apth_out), AcousticAnalogies.pressure_monopole(apth_out); marker=:x, label="AcousticAnalogies.combine")
scatter!(ax2_2, AcousticMetrics.time(apth_out), AcousticAnalogies.pressure_dipole(apth_out); marker=:x)
lines!(ax2_1, AcousticMetrics.time(apth_test), AcousticAnalogies.pressure_monopole(apth_test); label="Exact")
lines!(ax2_2, AcousticMetrics.time(apth_test), AcousticAnalogies.pressure_dipole(apth_test))
hidexdecorations!(ax2_1, grid=false)
axislegend(ax2_1; merge=true, unique=true, framevisible=false, bgcolor=:transparent, position=:rt)
save("combine_test.png", fig2)
nothing # hide
```
![](combine_test.png)

Right on top of each other.

### F1A Tests
The most complicated part of AcousticAnalogies.jl is the implementation of the F1A calculation itself.
For example, the compact form of the F1A dipole term as implemented in AcousticAnalogies.jl (neglecting surface deformation) is

```math
4 \pi c_0 p_d = \int_{L=0} \left[ \left( \dot{\vec{f}} \cdot \vec{D}_{1A} + \vec{f} \cdot \vec{E}_{1A} \right) dr \right]
```

where 

  * ``p_d`` is the "dipole" part of the acoustic pressure
  * ``\vec{f}`` is the loading per unit span on the source element
  * ``\vec{\dot{f}}`` is the source-time derivative of the loading per unit span on the source element
  * ``dr`` is the differential length of the source element
  * ``c_0`` is the ambient speed of sound
  * ``\vec{D}_{1A}`` and ``\vec{E}_{1A}`` are complicated functions of the position, velocity, and acceleration of the source element
  * ``L = 0`` indicates the integration is performed over some curve defined by ``L = 0``.

How are we going to test that we have all that implemented properly?
Well, it turns out that Farassat's original formulation (F1) is much simpler than F1A:

```math
4 \pi c_0 p_d = \frac{\partial}{\partial t} \int_{L=0} \left( \vec{f} \cdot \vec{B}_{1}\right) dr + \int_{L=0}\left( \vec{f} \cdot \vec{C}_1 \right) dr
```

where ``\vec{B}_1`` and ``\vec{C}_1`` are again functions of the position of the source element and time derivatives of the same.
It might not look that much simpler, but it is, because:

  * The F1 integrands don't depend on ``\dot{\vec{f}}``
  * ``\vec{D}_{1A}`` and ``\vec{E}_{1A}`` from F1A are more complicated than ``\vec{B}_1`` and ``\vec{C}_1``, and involve higher-order time derivatives

But the key thing to understand about F1 and F1A is that they are equivalent—going from F1A to F1 involves some fancy math (moving the derivative with respect to the observer time $t$ inside the integral), but should give the same answer.
The only trick is this: how will we evaluate the derivatives with respect to the observer time ``t`` in the F1 expressions?
What we'll do here is just use standard second-order-accurate finite difference approximations, i.e.,

```math
\frac{\partial g}{\partial t} = \frac{g(t+\Delta t) - g(t-\Delta t)}{2 \Delta t} + \mathcal{O}(\Delta t^2)
```

where the notation ``\mathcal{O}(\Delta t^2)`` indicates that the error associated with the finite difference approximation should be proportional to ``\Delta t^2``.
But this means that we can't expect our F1 calculation to exactly match F1A.
So, what to do about that?
What we can expect is that, if F1A and F1 have been implemented properly, the difference between them should go to zero at a second-order rate.
So we can systematically reduce the time step size ``\Delta t`` used to evaluate F1, and check that goes to zero at the expected rate.
If it does, that proves that the only error between the two codes is due to the finite difference approximation, and gives us strong evidence that both F1A and F1 have been implemented properly.

So, let's try it out!
First we'll need a function that evaluates the f1 integrands

```@example f1a_tests
using AcousticAnalogies
using LinearAlgebra: norm
using NLsolve
using Polynomials
using GLMakie

function f1_integrand(se, obs, t)
    c0 = se.c0

    # Need to get the retarded time.
    R(τ) = [t - (τ[1] + norm(obs(t) .- se.y0dot(τ[1]))/c0)]
    result = nlsolve(R, [-0.1], autodiff=:forward)
    if !converged(result)
        @error "nlsolve retarded time calculation did not converge:\n$(result)"
    end
    τ = result.zero[1]

    # Position of source at the retarted time.
    y = se.y0dot(τ)

    # Position vector from source to observer.
    rv = obs(t) .- y

    # Distance from source to observer.
    r = AcousticAnalogies.norm_cs_safe(rv)

    # Unit vector pointing from source to observer.
    rhat = rv./r

    # First time derivative of rv.
    rv1dot = -se.y1dot(τ)

    # Mach number of the velocity of the source in the direction of the
    # observer.
    Mr = AcousticAnalogies.dot_cs_safe(-rv1dot/se.c0, rhat)

    # Now evaluate the integrand.
    p_m_integrand = se.ρ0/(4*pi)*se.Λ*se.Δr/(r*(1 - Mr))

    # Loading at the retarded time.
    f0dot = se.f0dot(τ)

    p_d_integrand_ff = (1/(4*pi*c0))*AcousticAnalogies.dot_cs_safe(f0dot, rhat)/(r*(1 - Mr))*se.Δr
    p_d_integrand_nf = (1/(4*pi*c0))*AcousticAnalogies.dot_cs_safe(f0dot, rhat)*c0/(r^2*(1 - Mr))*se.Δr

    return τ, p_m_integrand, p_d_integrand_ff, p_d_integrand_nf
end

```

The `f1_integrand` function takes a source element `se`, and acoustic observer `obs`, and an observer time `t` and finds the source time and intermediate stuff that will eventually want to differentiate using the finite difference approximation.

Now we need to make up a source and observer that we can test this out with:
```@example f1a_tests

# https://docs.makie.org/v0.19/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

function doit()
    # Scale up the density to make the error bigger.
    rho = 1.226e6  # kg/m^3
    c0 = 340.0  # m/s
    Rtip = 1.1684  # meters
    radii = 0.99932*Rtip
    dradii = (0.99932 - 0.99660)*Rtip  # m
    area_over_chord_squared = 0.064
    chord = 0.47397E-02 * Rtip
    Λ = area_over_chord_squared * chord^2

    theta = 90.0*pi/180.0
    x0 = [cos(theta), 0.0, sin(theta)].*100.0.*12.0.*0.0254  # 100 ft in meters
    obs = StationaryAcousticObserver(x0)

    # Need the position and velocity of the source as a function of
    # source/retarded time. How do I want it to move? I want it to rotate around
    # an axis on the origin, pointing in the x direction.
    rpm = 2200
    omega = 2*pi/60*rpm
    period = 60/rpm
    fn = 180.66763939805125
    fc = 19.358679206883078
    y0dot(τ) = [0,  radii*cos(omega*τ), radii*sin(omega*τ)]
    y1dot(τ) = [0, -omega*radii*sin(omega*τ), omega*radii*cos(omega*τ)]
    y2dot(τ) = [0, -omega^2*radii*cos(omega*τ), -omega^2*radii*sin(omega*τ)]
    y3dot(τ) = [0, omega^3*radii*sin(omega*τ), -omega^3*radii*cos(omega*τ)]
    f0dot(τ) = [-fn, -sin(omega*τ)*fc, cos(omega*τ)*fc]
    f1dot(τ) = [0, -omega*cos(omega*τ)*fc, -omega*sin(omega*τ)*fc]
    u(τ) = y0dot(τ)./radii
    sef1 = CompactSourceElement(rho, c0, dradii, Λ, y0dot, y1dot, nothing, nothing, f0dot, nothing, 0.0, u)

    t = 0.0
    dt = period*0.5^4

    τ0, pmi0, pdiff0, pdinf0 = f1_integrand(sef1, obs, t)
    sef1a = CompactSourceElement(rho, c0, dradii, Λ, y0dot(τ0), y1dot(τ0), y2dot(τ0), y3dot(τ0), f0dot(τ0), f1dot(τ0), τ0, u(τ0))
    apth = f1a(sef1a, obs)

    err_prev_pm = nothing
    err_prev_pd = nothing
    dt_prev = nothing
    dt_curr = dt
    first_time = true

    err_pm = Vector{Float64}()
    err_pd = Vector{Float64}()
    dts = Vector{Float64}()
    ooa_pm = Vector{Float64}()
    ooa_pd = Vector{Float64}()
    # Gradually reduce time step size, recording the error and order-of-accuracy each time.
    for n in 1:7
        τ_1, pmi_1, pdiff_1, pdinf_1 = f1_integrand(sef1, obs, t-dt_curr)
        τ1, pmi1, pdiff1, pdinf1 = f1_integrand(sef1, obs, t+dt_curr)

        p_m_f1 = (pmi_1 - 2*pmi0 + pmi1)/(dt_curr^2)
        p_d_f1 = (pdiff1 - pdiff_1)/(2*dt_curr) + pdinf0

        err_curr_pm = abs(p_m_f1 - apth.p_m)
        err_curr_pd = abs(p_d_f1 - apth.p_d)

        if first_time
            first_time = false
        else
            push!(ooa_pm, log(err_curr_pm/err_prev_pm)/log(dt_curr/dt_prev))
            push!(ooa_pd, log(err_curr_pd/err_prev_pd)/log(dt_curr/dt_prev))
        end

        push!(dts, dt_curr)
        push!(err_pm, err_curr_pm)
        push!(err_pd, err_curr_pd)

        dt_prev = dt_curr
        err_prev_pm = err_curr_pm
        err_prev_pd = err_curr_pd
        dt_curr = 0.5*dt_curr
    end

    # Fit a line through the errors on a log-log plot, then check that the slope
    # is second-order.
    l = fit(log.(dts), log.(err_pm), 1)
    println("monopole term convergence rate = $(l.coeffs[2])")

    l = fit(log.(dts), log.(err_pd), 1)
    println("dipole term convergence rate = $(l.coeffs[2])")

    # Plot the error and observered order of accuracy.
    fig = Figure()
    ax1 = fig[1, 1] = Axis(fig, xlabel="time step size", ylabel="error", xscale=log10, xticks=LogTicks(IntegerTicks()), yscale=log10)
    ax2 = fig[2, 1] = Axis(fig, xlabel="time step size", ylabel="convergence rate", xscale=log10, xticks=LogTicks(IntegerTicks()))
    linkxaxes!(ax2, ax1)
    lines!(ax1, dts, err_pm, label="monopole term")
    lines!(ax1, dts, err_pd, label="dipole term")
    lines!(ax2, dts[2:end], ooa_pm, label="monopole term")
    lines!(ax2, dts[2:end], ooa_pd, label="dipole term")
    ylims!(ax2, -0.1, 3.1)
    axislegend(ax1; merge=true, unique=true, framevisible=false, bgcolor=:transparent, position=:lt)
    save("f1a_test.png", fig)
end

doit()
```
![](f1a_test.png)

The convergence rate of the error (the bottom plot) is extremely close to 2, which is what we're looking for.

### ANOPP2 Comparisons
The AcousticAnalogies.jl test suite includes comparisons to ANOPP2 predictions.
Tests for a hypothetical isolated rotor are performed over a range of RPMs, with both stationary and moving observers.
Here is an example using a moving observer:

```@example anopp2
using GLMakie
using FLOWMath: akima
include(joinpath(@__DIR__, "..", "..", "test", "anopp2_run.jl"))
using .ANOPP2Run
rpm = 2200.0
t, p_thickness, p_loading, p_monopole_a2, p_dipole_a2 = ANOPP2Run.get_results(;
    stationary_observer=false, theta=0.0, f_interp=akima, rpm=rpm, irpm=11)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig, xlabel="time, blade passes", ylabel="acoustic pressure, monopole, Pa")
ax2 = fig[2, 1] = Axis(fig, xlabel="time, blade passes", ylabel="acoustic pressure, dipole, Pa")
lines!(ax1, t, p_thickness, label="AcousticAnalogies.jl")
scatter!(ax1, t, p_monopole_a2, label="ANOPP2", markersize=6)
lines!(ax2, t, p_loading)
scatter!(ax2, t, p_dipole_a2, markersize=6)
hidexdecorations!(ax1, grid=false)
axislegend(ax1; merge=true, unique=true, framevisible=false, bgcolor=:transparent, position=:lt)
save("anopp2_comparison.png", fig)
```
![](anopp2_comparison.png)

The difference between the two codes' predictions is very small (less than 1% error).

### Brooks, Pope & Marcolini Tests
```@example bpm_bl_thickness
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/10^6", ylabel="δ_0/c",
                       xscale=log10, yscale=log10,
                       xminorticksvisible=true, yminorticksvisible=true,
                       xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()), yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure06-bl_thickness-untripped.csv")
bpm_untripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_untripped[:, 1]
deltastar0_c = bpm_untripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="untripped, BPM report", color=colors[2])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.bl_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="untripped, Julia", color=colors[2])

xlims!(ax1, 0.04, 3)
ylims!(ax1, 0.01, 0.2)
axislegend(ax1)
save("19890016302-figure06-bl_thickness.png", fig)
```
![](19890016302-figure06-bl_thickness.png)

```@example bpm_disp_thickness
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/10^6", ylabel="δ_0^*/c",
                       xscale=log10, yscale=log10,
                       xminorticksvisible=true, yminorticksvisible=true,
                       xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()), yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure06-disp_thickness-tripped.csv")
bpm_tripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_tripped[:, 1]
deltastar0_c = bpm_tripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="tripped, BPM report", color=colors[1])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="tripped, Julia", color=colors[1])

# deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)
# lines!(ax1, Re_c_1e6_jl, deltastar0_c_interp; linestyle=:dash, label="tripped, interp", color=colors[1])

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure06-disp_thickness-untripped.csv")
bpm_untripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_untripped[:, 1]
deltastar0_c = bpm_untripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="untripped, BPM report", color=colors[2])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="untripped, Julia", color=colors[2])

# deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)
# lines!(ax1, Re_c_1e6_jl, deltastar0_c_interp; linestyle=:dash, label="untripped, interp", color=colors[2])

xlims!(ax1, 0.04, 3)
ylims!(ax1, 0.001, 0.03)
axislegend(ax1)
save("19890016302-figure06-disp_thickness.png", fig)
```
![](19890016302-figure06-disp_thickness.png)

```@example bpm_disp_thickness_star_tripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ^*/δ_0^*",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure07-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies.disp_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure07-suction_side.csv")
bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_suction_side[:, 1]
deltastar_bpm = bpm_suction_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, markersize=4, color=colors[2], label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies.disp_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[2], label="suction side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[2], linestyle=:dash, label="suction side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("suction side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 200)
axislegend(ax1, position=:lt)
save("19890016302-figure07.png", fig)
```
![](19890016302-figure07.png)

```@example bpm_bl_thickness_star_untripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ/δ_0",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure08-bl_thickness-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies.bl_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 40)
axislegend(ax1, position=:lt)
save("19890016302-figure08-bl_thickness.png", fig)
```
![](19890016302-figure08-bl_thickness.png)

```@example bpm_disp_thickness_star_untripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ^*/δ_0^*",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure08-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies.disp_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure08-suction_side.csv")
bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_suction_side[:, 1]
deltastar_bpm = bpm_suction_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, markersize=4, color=colors[2], label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies.disp_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[2], label="suction side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[2], linestyle=:dash, label="suction side, interp")

# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("suction side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 200)
axislegend(ax1, position=:lt)
save("19890016302-figure08.png", fig)
```
![](19890016302-figure08.png)

```@example bpm_K_1
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c", ylabel="Peak scaled SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure77.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
Re_c_bpm = bpm[:, 1]
K_1_bpm = bpm[:, 2]
scatter!(ax1, Re_c_bpm, K_1_bpm, color=colors[1], markersize=8, label="BPM report")

Re_c_jl = range(minimum(Re_c_bpm), maximum(Re_c_bpm); length=50)
K_1_jl = AcousticAnalogies.K_1.(Re_c_jl)
lines!(ax1, Re_c_jl, K_1_jl, color=colors[1], label="Julia")

# Interpolate:
# K_1_interp = linear(Re_c_bpm, K_1_bpm, Re_c_jl)
# lines!(ax1, Re_c_jl, K_1_interp; color=colors[1], linestyle=:dash, label="interp")

# # Check error.
# vmin, vmax = extrema(K_1_bpm)
# err = abs.(K_1_jl .- K_1_interp)./(vmax - vmin)
# println("K_1_error =\n$(err)")

xlims!(ax1, 10^4, 10^7)
ylims!(ax1, 110.0, 150.0)
axislegend(ax1, position=:lt)
save("19890016302-figure77.png", fig)
```
![](19890016302-figure77.png)

```@example bpm_A
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal number ratio, St/St_peak", ylabel="Function A level, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure78-A_min.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
A = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, A, color=colors[1], markersize=8, label="A_min, BPM report")

St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
A_jl = AcousticAnalogies.A.(St_St_peak_jl, 9.5e4)
lines!(ax1, St_St_peak_jl, A_jl, color=colors[1], label="A_min, Julia")

# # Interpolate:
# A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, A_interp; color=colors[1], linestyle=:dash, label="A_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(A)
# err = abs.(A_jl .- A_interp)./(vmax - vmin)
# println("A_min error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure78-A_max.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
A = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, A, color=colors[2], markersize=8, label="A_max, BPM report")

St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
A_jl = AcousticAnalogies.A.(St_St_peak_jl, 8.58e5)
lines!(ax1, St_St_peak_jl, A_jl, color=colors[2], label="A_max, Julia")

# # Interpolate:
# A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, A_interp; color=colors[2], linestyle=:dash, label="A_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(A)
# err = abs.(A_jl .- A_interp)./(vmax - vmin)
# println("A_max error =\n$(err)")

xlims!(ax1, 0.1, 20)
ylims!(ax1, -20.0, 0.0)
axislegend(ax1, position=:lt)
save("19890016302-figure78-A.png", fig)
```
![](19890016302-figure78-A.png)

```@example bpm_B
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal number ratio, St/St_peak", ylabel="Function B level, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure78-B_min.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
B = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, B, color=colors[1], markersize=8, label="B_min, BPM report")

# St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
St_St_peak_jl = range(0.5, 2; length=50)
B_jl = AcousticAnalogies.B.(St_St_peak_jl, 9.5e4)
lines!(ax1, St_St_peak_jl, B_jl, color=colors[1], label="B_min, Julia")

# # Interpolate:
# B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, B_interp; color=colors[1], linestyle=:dash, label="B_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(B)
# err = abs.(B_jl .- B_interp)./(vmax - vmin)
# println("B_min error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure78-B_max.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
B = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, B, color=colors[2], markersize=8, label="B_max, BPM report")

# St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
St_St_peak_jl = range(0.2, 4; length=50)
B_jl = AcousticAnalogies.B.(St_St_peak_jl, 8.58e5)
lines!(ax1, St_St_peak_jl, B_jl, color=colors[2], label="B_max, Julia")

# # Interpolate:
# B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, B_interp; color=colors[2], linestyle=:dash, label="B_max, interp")
# 
# # Check error.
# vmin, vmax = extrema(B)
# err = abs.(B_jl .- B_interp)./(vmax - vmin)
# println("B_max error =\n$(err)")

xlims!(ax1, 0.1, 20)
ylims!(ax1, -20.0, 0.0)
axislegend(ax1, position=:lt)
save("19890016302-figure78-B.png", fig)
```
![](19890016302-figure78-B.png)

```@example bpm_St_2
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Angle of attack α^*, deg", ylabel="Peak Strouhal number, St_peak",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure80-M0.093.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
St_2 = bpm[:, 2]
scatter!(ax1, alpha_deg, St_2, color=colors[1], markersize=8, label="St_2 for M = 0.093, BPM")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.093), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, St_2_jl, color=colors[1], label="St_2 for M = 0.093, Julia")

# # Interpolate:
# St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, St_2_interp; color=colors[1], linestyle=:dash, label="St_2 for M = 0.093, interp")
# 
# # Check error.
# vmin, vmax = extrema(St_2)
# err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
# println("St_2 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure80-M0.209.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
St_2 = bpm[:, 2]
scatter!(ax1, alpha_deg, St_2, color=colors[2], markersize=8, label="St_2 for M = 0.209, BPM")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.209), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, St_2_jl, color=colors[2], label="St_2 for M = 0.209, Julia")

# # Interpolate:
# St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, St_2_interp; color=colors[2], linestyle=:dash, label="St_2 for M = 0.209, interp")
# 
# # Check error.
# vmin, vmax = extrema(St_2)
# err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
# println("St_2 error =\n$(err)")

xlims!(ax1, 0.0, 25.0)
ylims!(ax1, 0.01, 1)
axislegend(ax1, position=:lt)
save("19890016302-figure80.png", fig)
```
![](19890016302-figure80.png)

```@example bpm_K_2_K_1
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Angle of attack α_*, deg", ylabel="Extracted scaled levels minus K_1, dB")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure82-M0.093.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[1], markersize=8, label="M = 0.093, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.093, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[1], label="M = 0.093, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[1], linestyle=:dash, label="M = 0.093, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure82-M0.116.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[2], markersize=8, label="M = 0.116, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.116, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[2], label="M = 0.116, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[2], linestyle=:dash, label="M = 0.116, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure82-M0.163.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[3], markersize=8, label="M = 0.163, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.163, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[3], label="M = 0.163, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[3], linestyle=:dash, label="M = 0.163, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure82-M0.209.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[4], markersize=8, label="M = 0.209, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.209, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[4], label="M = 0.209, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[4], linestyle=:dash, label="M = 0.209, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

xlims!(ax1, 0.0, 25.0)
ylims!(ax1, -20, 20)
axislegend(ax1, position=:lt)
save("19890016302-figure82.png", fig)
```
![](19890016302-figure82.png)

```@example bpm_figure11_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure11-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0
alphastar0 = 12.5*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

# Now, I want to know what branches were actually used.
# @show getproperty.(tblte_branches, :K_1) getproperty.(tblte_branches, :DeltaK_1) getproperty.(tblte_branches, :A_s) getproperty.(tblte_branches, :A_p) getproperty.(tblte_branches, :B) getproperty.(tblte_branches, :St_2) getproperty.(tblte_branches, :K_2)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure11-a.png", fig)
```
![](19890016302-figure11-a.png)

```@example bpm_figure11_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure11-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
# M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
M = U/340.46
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0
# Re_c = U*chord/nu
# deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# deltastar_s = 0.0093*chord
# deltastar_s *= 0.6  # "lightly tripped" boundary layer?
# deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# deltastar_p *= 0.6  # "lightly tripped" boundary layer?
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
# St_s = @. f_jl*deltastar_s/U
# St_p = @. f_jl*deltastar_p/U
# # St_1 = AcousticAnalogies.St_1(M)
# St_1 = 0.056 # From Table 1 in the appendix.
# A_s = AcousticAnalogies.A.(St_s./St_1, Re_c)
# A_p = AcousticAnalogies.A.(St_p./St_1, Re_c)
# Re_deltastar_p = U*deltastar_p/nu
# ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
# # K_1 = AcousticAnalogies.K_1(Re_c)
# # K_1 = 128.0
# println("K_1 = $K_1 (should be 128.5)")
# SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
# SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
# @show M Re_c
# @show deltastar_s/chord
# @show AcousticAnalogies.St_1(M)

alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

# Now, I want to know what branches were actually used.
@show getproperty.(tblte_branches, :K_1) getproperty.(tblte_branches, :DeltaK_1) getproperty.(tblte_branches, :A_s) getproperty.(tblte_branches, :A_p) getproperty.(tblte_branches, :B) getproperty.(tblte_branches, :St_2) getproperty.(tblte_branches, :K_2)

# K_1_table = 128.0  # Value from line 11(b) in Table (1), BPM appendix.

# SPL_s_jl_table = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1_table - 3
# SPL_p_jl_table = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1_table - 3 + ΔK_1

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# lines!(ax1, f_jl./1e3, SPL_s_jl_table; label="TBL-TE suction side, Julia, K_1 from Table 1")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# lines!(ax1, f_jl./1e3, SPL_p_jl_table; label="TBL-TE pressure side, Julia, K_1 from Table 1")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure11-b.png", fig)
```
![](19890016302-figure11-b.png)

```@example bpm_figure11_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure11-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
# M = 0.116  # Mach number, corresponds to U = 36.6 m/s in BPM report
M = U/340.46
@show M
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0
Re_c = U*chord/nu
deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
St_s = @. f_jl*deltastar_s/U
St_p = @. f_jl*deltastar_p/U
A_s = AcousticAnalogies.A.(St_s./AcousticAnalogies.St_1(M), Re_c)
A_p = AcousticAnalogies.A.(St_p./AcousticAnalogies.St_1(M), Re_c)
Re_deltastar_p = U*deltastar_p/nu
ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
K_1 = AcousticAnalogies.K_1(Re_c)
SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
@show maximum(SPL_s_jl)
@show Re_c

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure11-c.png", fig)
```
![](19890016302-figure11-c.png)

```@example bpm_figure11_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure11-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0
Re_c = U*chord/nu
deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
St_s = @. f_jl*deltastar_s/U
St_p = @. f_jl*deltastar_p/U
A_s = AcousticAnalogies.A.(St_s./AcousticAnalogies.St_1(M), Re_c)
A_p = AcousticAnalogies.A.(St_p./AcousticAnalogies.St_1(M), Re_c)
Re_deltastar_p = U*deltastar_p/nu
ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
K_1 = AcousticAnalogies.K_1(Re_c)
SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
@show maximum(SPL_s_jl)
@show Re_c

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure11-d.png", fig)
```
![](19890016302-figure11-d.png)

```@example bpm_figure12_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-U71.3-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-U71.3-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-U71.3-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 1.5*pi/180
# D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
# Re_c = U*chord/nu
# deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# f = range(0.2, 20; length=100).*1e3
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
# St_s = @. f_jl*deltastar_s/U
# St_p = @. f_jl*deltastar_p/U
# A_s = AcousticAnalogies.A.(St_s./AcousticAnalogies.St_1(M), Re_c)
# A_p = AcousticAnalogies.A.(St_p./AcousticAnalogies.St_1(M), Re_c)
# Re_deltastar_p = U*deltastar_p/nu
# ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
# K_1 = AcousticAnalogies.K_1(Re_c)
# SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
alphastar0 = 12.5*pi/180
SPL_s_jl = AcousticAnalogies.TBL_TE_s.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
# SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
SPL_p_jl = AcousticAnalogies.TBL_TE_p.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
# B = AcousticAnalogies.B.(St_s./AcousticAnalogies.St_2(AcousticAnalogies.St_1(M), alphastar), Re_c)
# K_2 = AcousticAnalogies.K_2(Re_c, M, alphastar)
# SPL_alpha_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + B + K_2
SPL_alpha_jl = AcousticAnalogies.TBL_TE_alpha.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
gamma0 = 23.43*M + 4.651
println("gamma0 = $gamma0 deg, alphastar = $(alphastar*180/pi) deg")

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 12 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure12-a.png", fig)
```
![](19890016302-figure12-a.png)

```@example bpm_figure12_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-b-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure12-b-TBL-TE-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 36.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 1.5*pi/180
# D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
# Re_c = U*chord/nu
# deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.TrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
f_jl = ExactThirdOctaveCenterBands(0.2, 20e3)
# St_s = @. f_jl*deltastar_s/U
# St_p = @. f_jl*deltastar_p/U
# A_s = AcousticAnalogies.A.(St_s./AcousticAnalogies.St_1(M), Re_c)
# A_p = AcousticAnalogies.A.(St_p./AcousticAnalogies.St_1(M), Re_c)
# Re_deltastar_p = U*deltastar_p/nu
# ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
# K_1 = AcousticAnalogies.K_1(Re_c)
# SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
alphastar0 = 12.5*pi/180
SPL_s_jl = AcousticAnalogies.TBL_TE_s.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
# SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
SPL_p_jl = AcousticAnalogies.TBL_TE_p.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
SPL_alpha_jl = AcousticAnalogies.TBL_TE_alpha.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))
@show maximum(SPL_s_jl)
# @show Re_c

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 12 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure12-b.png", fig)
```
![](19890016302-figure12-b.png)

```@example bpm_figure26_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure26-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure26-a.png", fig)
```
![](19890016302-figure26-a.png)

```@example bpm_figure26_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure26-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure26-b.png", fig)
```
![](19890016302-figure26-b.png)

```@example bpm_figure26_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure26-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure26-c.png", fig)
```
![](19890016302-figure26-c.png)

```@example bpm_figure26_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure26-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure26-d.png", fig)
```
![](19890016302-figure26-d.png)

```@example bpm_figure28_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-a-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure28-a.png", fig)
```
![](19890016302-figure28-a.png)

```@example bpm_figure28_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-b-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-b-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-b.png", fig)
```
![](19890016302-figure28-b.png)

```@example bpm_figure28_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-c-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-c-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-c.png", fig)
```
![](19890016302-figure28-c.png)

```@example bpm_figure28_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-d-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure28-d-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, khz", ylabel="spl_1/3, db",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="figure 28 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, bpm")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, bpm")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, bpm")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-d.png", fig)
```
![](19890016302-figure28-d.png)

```@example bpm_figure45_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure45-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure45-a-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure45-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure45-a-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
# M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
M = U/340.46
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
# D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 1.5*pi/180
# Re_c = U*chord/nu
# deltastar_s = AcousticAnalogies.disp_thickness_s(AcousticAnalogies.UntrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
# deltastar_p = AcousticAnalogies.disp_thickness_p(AcousticAnalogies.UntrippedN0012BoundaryLayer(), Re_c, alphastar)*chord
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
# St_s = @. f_jl*deltastar_s/U
# St_p = @. f_jl*deltastar_p/U
# A_s = AcousticAnalogies.A.(St_s./AcousticAnalogies.St_1(M), Re_c)
# A_p = AcousticAnalogies.A.(St_p./AcousticAnalogies.St_1(M), Re_c)
# Re_deltastar_p = U*deltastar_p/nu
# ΔK_1 = AcousticAnalogies.DeltaK_1(alphastar, Re_deltastar_p)
# K_1 = AcousticAnalogies.K_1(Re_c)
# SPL_s_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + A_s + K_1 - 3
# SPL_p_jl = @. 10*log10((deltastar_p*M^5*L*D_h)/(r_e^2)) + A_p + K_1 - 3 + ΔK_1
# B = AcousticAnalogies.B.(St_s./AcousticAnalogies.St_2(AcousticAnalogies.St_1(M), alphastar), Re_c)
# K_2 = AcousticAnalogies.K_2(Re_c, M, alphastar)
# SPL_alpha_jl = @. 10*log10((deltastar_s*M^5*L*D_h)/(r_e^2)) + B + K_2
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 45 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker='◇', label="LBL-VS, BPM")
scatterlines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; marker='◇', label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure45-a.png", fig)
```
![](19890016302-figure45-a.png)

```@example bpm_figure48_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure48-c-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure48-c-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure48-c-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure48-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 22.86e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 48 (c) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure48-c.png", fig)
```
![](19890016302-figure48-c.png)

```@example bpm_figure54_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure54-a-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure54-a-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure54-a-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure54-a-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 15.24e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 2.7*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 54 (a) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 50, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure54-a.png", fig)
```
![](19890016302-figure54-a.png)

```@example bpm_figure59_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure59-c-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure59-c-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure59-c-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure59-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 59 (c) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure59-c.png", fig)
```
![](19890016302-figure59-c.png)

```@example bpm_figure60_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-c-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-c-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-c-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 3.3*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 60 (c) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure60-c.png", fig)
```
![](19890016302-figure60-c.png)

```@example bpm_figure60_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-d-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-d-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-d-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure60-d-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 3.3*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 60 (d) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure60-d.png", fig)
```
![](19890016302-figure60-d.png)

```@example bpm_figure65_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure65-d-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure65-d-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure65-d-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure65-d-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 65 (d) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 50, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure65-d.png", fig)
```
![](19890016302-figure65-d.png)

```@example bpm_figure66_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure66-b-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure66-b-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure66-b-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure66-b-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 4.2*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_lbl_vs_jl = AcousticAnalogies.LBL_VS.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 66 (b) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
# lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")
# 
# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
scatterlines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure66-b.png", fig)
```
![](19890016302-figure66-b.png)

```@example bpm_figure69_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-a-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-a-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 15.4*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 69 (a) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 60, 100)
axislegend(ax1, position=:rt)
save("19890016302-figure69-a.png", fig)
```
![](19890016302-figure69-a.png)

```@example bpm_figure69_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-b-TBL-TE-suction.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_s = bpm[:, 1]
# SPL_s = bpm[:, 2]

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-b-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure69-b-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 15.4*pi/180
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
alphastar0 = 12.5*pi/180
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)
SPL_alpha_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 3)
tblte_branches = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 4)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 69 (b) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure69-b.png", fig)
```
![](19890016302-figure69-b.png)

```@example bpm_St_1_prime
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c", ylabel="Peak Strouhal number, St'_peak",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure86-St_1_prime.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
Re_c_bpm = bpm[:, 1]
St_1_prime_bpm = bpm[:, 2]
scatter!(ax1, Re_c_bpm, St_1_prime_bpm, color=colors[1], markersize=4, label="BPM")

Re_c_jl = 10.0.^(range(4, 7; length=100))
St_1_prime_jl = AcousticAnalogies.St_1_prime.(Re_c_jl)
lines!(ax1, Re_c_jl, St_1_prime_jl, color=colors[1], label="Julia")

xlims!(ax1, 1e4, 1e7)
ylims!(ax1, 0.01, 1)
axislegend(ax1, position=:lt)
save("19890016302-figure86.png", fig)
```
![](19890016302-figure86.png)

```@example bpm_lbl_vs_G1
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="St'/St'_peak", ylabel="Function G_1 level, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure85-G1.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
e_bpm = bpm[:, 1]
G1_bpm = bpm[:, 2]
scatter!(ax1, e_bpm, G1_bpm, color=colors[1], markersize=4, label="BPM")

e_jl = 10.0.^(range(-1, 1; length=101))
G1_jl = AcousticAnalogies.G1.(e_jl)

lines!(ax1, e_jl, G1_jl, colors=colors[1], label="Julia")

xlims!(ax1, 0.1, 10)
ylims!(ax1, -30, 0)
axislegend(ax1, position=:lt)
save("19890016302-figure85.png", fig)
```
![](19890016302-figure85.png)

```@example bpm_lbl_vs_St_peak_prime_alphastar
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha^*, deg", ylabel="St'_peak/St'_1",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure87.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alphastar_bpm = bpm[:, 1]
St_peak_ratio_bpm = bpm[:, 2]
scatter!(ax1, alphastar_bpm, St_peak_ratio_bpm, color=colors[1], markersize=4, label="BPM")

St_1_prime = 0.25  # Just make up a value, since we're multiplying and then dividing by it anyway.
alphastar_jl = range(0.0*pi/180, 7.0*pi/180; length=21)
St_peak_ratio_jl = AcousticAnalogies.St_peak_prime.(St_1_prime, alphastar_jl)./St_1_prime
lines!(ax1, alphastar_jl.*180/pi, St_peak_ratio_jl, color=colors[1], label="Julia")

xlims!(ax1, 0, 7)
ylims!(ax1, 0.5, 2)
axislegend(ax1, position=:lt)
save("19890016302-figure87.png", fig)
```
![](19890016302-figure87.png)

```@example bpm_lbl_vs_G2_alphastar
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/Re_c0", ylabel="G2 + G3",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure88-G2-alpha0.csv")
alphastar = 0.0*pi/180
bpm = DelimitedFiles.readdlm(fname, ',')
Re_c_bpm = bpm[:, 1]
G2_bpm = bpm[:, 2]
scatter!(ax1, Re_c_bpm, G2_bpm, color=colors[1], markersize=4, label="BPM - α^* = $(alphastar*180/pi)°")

Re_c_jl = 10.0.^range(log10(first(Re_c_bpm)), log10(last(Re_c_bpm)), length=51)
Re_c0 = AcousticAnalogies.Re_c0(alphastar)
Re_ratio_jl = Re_c_jl./Re_c0
# G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ 171.04 .- 3.03*(alphastar*180/pi)
G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ AcousticAnalogies.G3.(alphastar)
lines!(ax1, Re_c_jl, G2_jl, color=colors[1], label="Julia - α^* = $(alphastar*180/pi)°")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure88-G2-alpha6.csv")
alphastar = 6.0*pi/180
bpm = DelimitedFiles.readdlm(fname, ',')
Re_c_bpm = bpm[:, 1]
G2_bpm = bpm[:, 2]
scatter!(ax1, Re_c_bpm, G2_bpm, color=colors[2], markersize=4, label="BPM - α^* = $(alphastar*180/pi)°")

Re_c_jl = 10.0.^range(log10(first(Re_c_bpm)), log10(last(Re_c_bpm)), length=51)
Re_c0 = AcousticAnalogies.Re_c0(alphastar)
Re_ratio_jl = Re_c_jl./Re_c0
# G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ 171.04 .- 3.03*(alphastar*180/pi)
G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ AcousticAnalogies.G3.(alphastar)
lines!(ax1, Re_c_jl, G2_jl, color=colors[2], label="Julia - α^* = $(alphastar*180/pi)°")

xlims!(ax1, 10^4, 10^7)
ylims!(ax1, 125, 175)
axislegend(ax1, position=:lt)
save("19890016302-figure88.png", fig)
```
![](19890016302-figure88.png)

```@example bpm_lbl_vs_G2
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/Re_c0", ylabel="G2",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure89.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
Re_ratio_bpm = bpm[:, 1]
G2_bpm = bpm[:, 2]
scatter!(ax1, Re_ratio_bpm, G2_bpm, color=colors[1], markersize=4, label="BPM")

Re_ratio_jl = 10.0.^range(-1, 1, length=51)
G2_jl = AcousticAnalogies.G2.(Re_ratio_jl)
lines!(ax1, Re_ratio_jl, G2_jl, color=colors[1], label="Julia")

xlims!(ax1, 0.1, 100)
ylims!(ax1, -45, 5)
axislegend(ax1, position=:lt)
save("19890016302-figure89.png", fig)
```
![](19890016302-figure89.png)

```@example bpm_figure91
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure91-tip.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_tip = bpm[:, 1]
SPL_tip = bpm[:, 2]

# nu = 1.4529e-5  # kinematic viscosity, m^2/s
# L = 30.48e-2  # span in meters
chord = 15.24e-2  # chord in meters
speedofsound = 340.46
U = 71.3  # freestream velocity in m/s
# M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
M = U/speedofsound
M_c = 0.8*M
# speedofsound = U/M
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphatip = 10.8*pi/180
alphatip_prime = 0.71*alphatip
# Equation 64 in the BPM report.
M_max = (1 + 0.036*(alphatip_prime*180/pi))*M
U_max = M_max*speedofsound
f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_tip_jl = AcousticAnalogies.TIP.(f_jl, chord, M, M_c, U_max, M_max, r_e, θ_e, Φ_e, alphatip_prime, Ref(AcousticAnalogies.RoundedTip()))

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 91")
scatter!(ax1, f_tip, SPL_tip; marker='o', label="Tip, BPM")
lines!(ax1, f_jl./1e3, SPL_tip_jl; label="Tip, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure91.png", fig)
```
![](19890016302-figure91.png)


```@example bpm_figure95
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure95-0Psi.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
h_over_deltastar_0Psi = bpm[:, 1]
St_3prime_peak_0Psi = bpm[:, 2]
 
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure95-14Psi.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
h_over_deltastar_14Psi = bpm[:, 1]
St_3prime_peak_14Psi = bpm[:, 2]

h_over_deltastar_jl = 10.0.^(range(-1, 1; length=51))
St_3prime_peak_0Psi_jl = AcousticAnalogies.St_3prime_peak.(h_over_deltastar_jl, 0.0*pi/180)
St_3prime_peak_14Psi_jl = AcousticAnalogies.St_3prime_peak.(h_over_deltastar_jl, 14.0*pi/180)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Thickness ratio, h/δ^*", ylabel="Peak Strouhal number, St'''_peak",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks()),
                       title="Figure 95")

scatter!(ax1, h_over_deltastar_0Psi, St_3prime_peak_0Psi; marker='o', label="Ψ = 0°, BPM")
lines!(ax1, h_over_deltastar_jl, St_3prime_peak_0Psi_jl; label="Ψ = 0°, Julia")

scatter!(ax1, h_over_deltastar_14Psi, St_3prime_peak_14Psi; marker='o', label="Ψ = 14°, BPM")
lines!(ax1, h_over_deltastar_jl, St_3prime_peak_14Psi_jl; label="Ψ = 14°, Julia")

xlims!(ax1, 0.2, 10.0)
ylims!(ax1, 0.05, 0.3)
axislegend(ax1, position=:rt)
save("19890016302-figure95.png", fig)
```
![](19890016302-figure95.png)

```@example bpm_figure96
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure96-0Psi.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
h_over_deltastar_0Psi = bpm[:, 1]
G4_0Psi = bpm[:, 2]
 
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure96-14Psi.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
h_over_deltastar_14Psi = bpm[:, 1]
G4_14Psi = bpm[:, 2]

h_over_deltastar_jl = 10.0.^(range(-1, 1; length=51))
G4_0Psi_jl = AcousticAnalogies.G4.(h_over_deltastar_jl, 0.0*pi/180)
G4_14Psi_jl = AcousticAnalogies.G4.(h_over_deltastar_jl, 14.0*pi/180)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Thickness ratio, h/δ^*", ylabel="Scaled peak SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 96")

scatter!(ax1, h_over_deltastar_0Psi, G4_0Psi; marker='o', label="Ψ = 0°, BPM")
lines!(ax1, h_over_deltastar_jl, G4_0Psi_jl; label="Ψ = 0°, Julia")

scatter!(ax1, h_over_deltastar_14Psi, G4_14Psi; marker='o', label="Ψ = 14°, BPM")
lines!(ax1, h_over_deltastar_jl, G4_14Psi_jl; label="Ψ = 14°, Julia")

xlims!(ax1, 0.1, 10.0)
ylims!(ax1, 110, 180)
axislegend(ax1, position=:lt)
save("19890016302-figure96.png", fig)
```
![](19890016302-figure96.png)

```@example bpm_figure97a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar0p25.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p25 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg0p25 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar0p43.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p43 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg0p43 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar0p50.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p50 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg0p50 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar0p54.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p54 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg0p54 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar0p62.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p62 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg0p62 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi14-h_over_deltastar1p20.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_1p20 = bpm[:, 1]
G5_14Psi_h_over_deltastar_avg1p20 = bpm[:, 2]

St_3prime_over_St_3prime_peak_jl = 10.0.^(range(-1, 10; length=1001))
G5_14Psi_h_over_deltastar_avg0p25_jl = AcousticAnalogies.G5_Psi14.(0.25, St_3prime_over_St_3prime_peak_jl)
G5_14Psi_h_over_deltastar_avg0p43_jl = AcousticAnalogies.G5_Psi14.(0.43, St_3prime_over_St_3prime_peak_jl)
G5_14Psi_h_over_deltastar_avg0p50_jl = AcousticAnalogies.G5_Psi14.(0.50, St_3prime_over_St_3prime_peak_jl)
G5_14Psi_h_over_deltastar_avg0p54_jl = AcousticAnalogies.G5_Psi14.(0.54, St_3prime_over_St_3prime_peak_jl)
G5_14Psi_h_over_deltastar_avg0p62_jl = AcousticAnalogies.G5_Psi14.(0.62, St_3prime_over_St_3prime_peak_jl)
G5_14Psi_h_over_deltastar_avg1p20_jl = AcousticAnalogies.G5_Psi14.(1.20, St_3prime_over_St_3prime_peak_jl)
 
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal ratio, St'''/St'''_peak", ylabel="G_5, Ψ=14°",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 97a")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p25, G5_14Psi_h_over_deltastar_avg0p25; label="h/δ^* = 0.25, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p25_jl; label="h/δ^* = 0.25, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p43, G5_14Psi_h_over_deltastar_avg0p43; label="h/δ^* = 0.43, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p43_jl; label="h/δ^* = 0.43, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p50, G5_14Psi_h_over_deltastar_avg0p50; label="h/δ^* = 0.50, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p50_jl; label="h/δ^* = 0.50, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p54, G5_14Psi_h_over_deltastar_avg0p54; label="h/δ^* = 0.54, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p54_jl; label="h/δ^* = 0.54, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p62, G5_14Psi_h_over_deltastar_avg0p62; label="h/δ^* = 0.62, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p62_jl; label="h/δ^* = 0.62, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_1p20, G5_14Psi_h_over_deltastar_avg1p20; label="h/δ^* = 1.20, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg1p20_jl; label="h/δ^* = 1.20, Julia")

xlims!(ax1, 0.1, 10.0)
ylims!(ax1, -30, 10)
axislegend(ax1, position=:rt)
save("19890016302-figure97a.png", fig)
```
![](19890016302-figure97a.png)

```@example bpm_figure97b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar0p25.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p25 = bpm[:, 1]
G5_0Psi_h_over_deltastar_avg0p25 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar0p43.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p43 = bpm[:, 1]
G5_0Psi_h_over_deltastar_avg0p43 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar0p50.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p50 = bpm[:, 1]
G5_0Psi_h_over_deltastar_avg0p50 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar0p54.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_0p54 = bpm[:, 1]
G5_0Psi_h_over_deltastar_avg0p54 = bpm[:, 2]

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar0p62.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# St_3prime_over_St_3prime_peak_0p62 = bpm[:, 1]
# G5_0Psi_h_over_deltastar_avg0p62 = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure97-Psi0-h_over_deltastar1p20.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_3prime_over_St_3prime_peak_1p20 = bpm[:, 1]
G5_0Psi_h_over_deltastar_avg1p20 = bpm[:, 2]

St_3prime_over_St_3prime_peak_jl = 10.0.^(range(-1, 10; length=1001))
G5_0Psi_h_over_deltastar_avg0p25_jl = AcousticAnalogies.G5_Psi0.(0.25, St_3prime_over_St_3prime_peak_jl)
G5_0Psi_h_over_deltastar_avg0p43_jl = AcousticAnalogies.G5_Psi0.(0.43, St_3prime_over_St_3prime_peak_jl)
G5_0Psi_h_over_deltastar_avg0p50_jl = AcousticAnalogies.G5_Psi0.(0.50, St_3prime_over_St_3prime_peak_jl)
G5_0Psi_h_over_deltastar_avg0p54_jl = AcousticAnalogies.G5_Psi0.(0.54, St_3prime_over_St_3prime_peak_jl)
# G5_0Psi_h_over_deltastar_avg0p62_jl = AcousticAnalogies.G5_Psi0.(0.62, St_3prime_over_St_3prime_peak_jl)
G5_0Psi_h_over_deltastar_avg1p20_jl = AcousticAnalogies.G5_Psi0.(1.20, St_3prime_over_St_3prime_peak_jl)
 
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal ratio, St'''/St'''_peak", ylabel="G_5, Ψ=0°",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 97b")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p25, G5_0Psi_h_over_deltastar_avg0p25; label="h/δ^* = 0.25, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p25_jl; label="h/δ^* = 0.25, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p43, G5_0Psi_h_over_deltastar_avg0p43; label="h/δ^* = 0.43, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p43_jl; label="h/δ^* = 0.43, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p50, G5_0Psi_h_over_deltastar_avg0p50; label="h/δ^* = 0.50, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p50_jl; label="h/δ^* = 0.50, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_0p54, G5_0Psi_h_over_deltastar_avg0p54; label="h/δ^* = 0.54, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p54_jl; label="h/δ^* = 0.54, Julia")

# scatter!(ax1, St_3prime_over_St_3prime_peak_0p62, G5_0Psi_h_over_deltastar_avg0p62; label="h/δ^* = 0.62, BPM", marker='o')
# lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p62_jl; label="h/δ^* = 0.62, Julia")

scatter!(ax1, St_3prime_over_St_3prime_peak_1p20, G5_0Psi_h_over_deltastar_avg1p20; label="h/δ^* = 1.20, BPM", marker='o')
lines!(ax1, St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg1p20_jl; label="h/δ^* = 1.20, Julia")

xlims!(ax1, 0.1, 10.0)
ylims!(ax1, -30, 10)
axislegend(ax1, position=:rt)
save("19890016302-figure97b.png", fig)
```
![](19890016302-figure97b.png)

```@example bpm_figure98_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-b-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 1.1e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure98-b.png", fig)
```
![](19890016302-figure98-b.png)

```@example bpm_figure98_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-c-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 1.9e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure98-c.png", fig)
```
![](19890016302-figure98-c.png)

```@example bpm_figure98_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure98-d-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 2.5e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure98-d.png", fig)
```
![](19890016302-figure98-d.png)

```@example bpm_figure99_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 1.1e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-b.png", fig)
```
![](19890016302-figure99-b.png)

```@example bpm_figure99_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-c-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 1.9e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-c.png", fig)
```
![](19890016302-figure99-c.png)

```@example bpm_figure99_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "19890016302-figure99-d-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 2.5e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180
alphastar0 = 12.5*pi/180

f_jl = ExactThirdOctaveCenterBands(0.2e3, 20e3)
SPL_s_SPL_p_SPL_alpha_branches = AcousticAnalogies.TBL_TE_branch.(f_jl, nu, L, chord, U, M, M_c, r_e, θ_e, Φ_e, alphastar, alphastar0, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))

SPL_s_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 1)
SPL_p_jl = getindex.(SPL_s_SPL_p_SPL_alpha_branches, 2)

SPL_teb_vs_jl = AcousticAnalogies.BLUNT.(f_jl, nu, L, chord, h, Psi, U, M, M_c, r_e, θ_e, Φ_e, alphastar, Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()))


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-d.png", fig)
```
![](19890016302-figure99-d.png)

## Signed Commits
The AcousticAnalogies.jl GitHub repository requires all commits to the `main` branch to be signed.
See the [GitHub docs on signing commits](https://docs.github.com/en/authentication/managing-commit-signature-verification/about-commit-signature-verification) for more information.

## Reporting Bugs
Users can use the [GitHub Issues](https://docs.github.com/en/issues/tracking-your-work-with-issues/about-issues) feature to report bugs and submit feature requests.
