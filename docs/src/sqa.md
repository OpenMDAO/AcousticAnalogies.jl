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

But the key thing to understand about F1 and F1A is that they are equivalent—going from F1A to F1 involves some fancy math (moving the derivative with respect to the observer time $t$ into the integral), but should give the same answer.

## Signed Commits

## Reporting Bugs
