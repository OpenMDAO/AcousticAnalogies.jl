using AcousticAnalogies
using BenchmarkTools
using DelimitedFiles
using FLOWMath: linear, akima
using Interpolations
using KinematicCoordinateTransformations
using LinearAlgebra: ×
using Printf: @sprintf
using StaticArrays

include(joinpath(@__DIR__, "../test/gen_test_data/gen_ccblade_data/constants.jl"))

const paramsfile = joinpath(@__DIR__, "params-current.json")
const resultsfile = joinpath(@__DIR__, "results-current.json")

function get_dradii(radii, Rhub, Rtip)
    # How do I get the radial spacing? Well, for the inner elements, I'll just
    # assume that the interfaces are midway in between the centers.
    r_interface = 0.5.*(radii[1:end-1] .+ radii[2:end])
    # Then just say that the blade begins at Rhub, and ends at Rtip.
    r_interface = vcat([Rhub], r_interface, [Rtip])
    # And now the distance between interfaces is the spacing.
    dradii = r_interface[2:end] .- r_interface[1:end-1]
    return dradii
end

function run_current(; load_params=true, save_params=false)
    rpm = 2200.0  # rev/min
    omega = rpm*(2*pi/60.0)

    # Get the normal and circumferential loading from the CCBlade output.
    i = 11
    fname = joinpath(@__DIR__, "../test/gen_test_data/gen_ccblade_data/ccblade_omega$(@sprintf "%02d" i).csv")
    data = DelimitedFiles.readdlm(fname, ',')
    fn = data[:, 1]
    fc = data[:, 2]

    # Blade passing period.
    bpp = 2*pi/omega/num_blades

    num_src_times = 256
    num_obs_times = 2*num_src_times
    num_radial = length(radii)
    num_sources = num_blades*num_radial
    obs_time_range = 4.0*bpp

    dradii = get_dradii(radii, Rhub, Rtip)
    cs_area =  area_over_chord_squared .* chord.^2

    # Observer angle, rad. Zero is sideline (in the rotor plane of rotation).
    theta = 0.0
    x0 = [cos(theta)*100*12*0.0254, 0.0, sin(theta)*100*12*0.0254]  # 100 ft in meters

    # For the moving observer, is at x0 at time t0, moving with constant
    # velocity v0_hub.
    t0 = 0.0
    v0_hub = v.*[0.0, 0.0, 1.0]

    obs_stationary = StationaryAcousticObserver(SVector{3}(x0))
    obs_moving = ConstVelocityAcousticObserver(t0, SVector{3}(x0), SVector{3}(v0_hub))

    apth = Array{AcousticPressure{Float64, Float64, Float64}, 3}(undef, num_src_times, num_radial, num_blades)
    apth_total = AcousticPressure(zeros(Float64, num_obs_times), zeros(Float64, num_obs_times), zeros(Float64, num_obs_times))

    linear_interpolations_jl(t_cp, p_cp, t) = interpolate((t_cp,), p_cp, Gridded(Linear())).(t)

    suite = BenchmarkGroup()

    s_s = suite["stationary"] = BenchmarkGroup()
    s_s["linear"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_stationary, $apth, $apth_total, $linear)
    s_s["akima"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_stationary, $apth, $apth_total, $akima)
    s_s["linear_interpolations_jl"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_stationary, $apth, $apth_total, $linear_interpolations_jl)

    s_m = suite["moving"] = BenchmarkGroup()
    s_m["linear"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_moving, $apth, $apth_total, $linear)
    s_m["akima"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_moving, $apth, $apth_total, $akima)
    s_m["linear_interpolations_jl"] = @benchmarkable run_cf1a($num_blades, $v, $omega, $radii, $dradii, $cs_area, $fn, $fc, $obs_time_range, $obs_moving, $apth, $apth_total, $linear_interpolations_jl)

    if load_params && isfile(paramsfile)
        # Load the benchmark parameters.
        # https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/doc/manual.md#caching-parameters
        loadparams!(suite, BenchmarkTools.load(paramsfile)[1])

        # Also need to warmup the benchmarks to get rid of the JIT overhead
        # (when not using tune!):
        # https://discourse.julialang.org/t/benchmarktools-theory-and-practice/5728
        warmup(suite, verbose=false)
    else
        tune!(suite, verbose=false)
    end

    results = run(suite, verbose=false)

    if save_params
        BenchmarkTools.save(paramsfile, params(suite))
    end

    return suite, results
end

function run_cf1a(num_blades, v, omega, radii, dradii, cs_area, fn, fc, obs_time_range, obs, apth, apth_total, f_interp)
    # This is the same as kinematic_trans_pipe, except the un-transformed SourceElement2
    # objects are never assigned to an intermediate variable.
    t0 = 0.0
    rot_axis = @SVector [0.0, 0.0, 1.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    v0_hub = v.*rot_axis

    num_radial = length(radii)
    num_src_times = size(apth, 1)

    # Blade Passing Period.
    bpp = 2*pi/omega/num_blades
    src_time_range = 5.0*bpp

    rot_trans = SteadyRotXTransformation(t0, omega, 0.0)
    global_trans = ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # This is just an array of the angular offsets of each blade.
    θs = 2*pi/num_blades.*(0:(num_blades-1))

    dt = src_time_range/(num_src_times - 1)
    src_times = t0 .+ (0:num_src_times-1).*dt

    # Reshape for broadcasting.
    θs = reshape(θs, 1, 1, num_blades)
    radii = reshape(radii, 1, num_radial, 1)
    dradii = reshape(dradii, 1, num_radial, 1)
    cs_area = reshape(cs_area, 1, num_radial, 1)
    fn = reshape(fn, 1, num_radial, 1)
    fc = reshape(fc, 1, num_radial, 1)
    src_times = reshape(src_times, num_src_times, 1, 1)  # This isn't really necessary.

    # Get all the transformations!
    trans = compose.(src_times, Ref(const_vel_trans), compose.(src_times, Ref(global_trans), Ref(rot_trans)))

    # Transform the source elements.
    ses = CompactSourceElement.(rho, c0, radii, θs, dradii, cs_area, fn, fc, src_times) .|> trans

    # Do the acoustics.
    apth .= f1a.(ses, Ref(obs))

    # Get the common observer time.
    common_obs_time!(apth_total.t, apth, obs_time_range, 1)

    # Combine all the sources into one acoustic pressure time history.
    combine!(apth_total, apth, 1; f_interp=f_interp)

    return apth_total.t, apth_total.p_m, apth_total.p_d
end

# function compare(; load_params=true, save_params=false, save_results=false)
function compare(; load_params=true, save_params=false, save_results=false)
    suite, results_new = run_current(load_params=load_params, save_params=save_params)

    results_old = BenchmarkTools.load(resultsfile)[1]

    println("Stationary observer, FLOWMath linear interpolation:")
    rold = results_old["stationary"]["linear"]
    rnew = results_new["stationary"]["linear"]
    display(judge(median(rnew), median(rold)))

    println("Stationary observer, FLOWMath Akima spline interpolation:")
    rold = results_old["stationary"]["akima"]
    rnew = results_new["stationary"]["akima"]
    display(judge(median(rnew), median(rold)))

    println("Stationary observer, Interpolations.jl linear interpolation:")
    rold = results_old["stationary"]["linear_interpolations_jl"]
    rnew = results_new["stationary"]["linear_interpolations_jl"]
    display(judge(median(rnew), median(rold)))

    println("Moving observer, FLOWMath linear interpolation:")
    rold = results_old["moving"]["linear"]
    rnew = results_new["moving"]["linear"]
    display(judge(median(rnew), median(rold)))

    println("Moving observer, FLOWMath Akima spline interpolation:")
    rold = results_old["moving"]["akima"]
    rnew = results_new["moving"]["akima"]
    display(judge(median(rnew), median(rold)))

    println("Moving observer, Interpolations.jl linear interpolation:")
    rold = results_old["moving"]["linear_interpolations_jl"]
    rnew = results_new["moving"]["linear_interpolations_jl"]
    display(judge(median(rnew), median(rold)))

    if save_results
        BenchmarkTools.save(resultsfile, results_new)
    end

    return suite, results_old, results_new
end

if !isinteractive()
    compare(; load_params=true, save_params=false, save_results=false)
end
