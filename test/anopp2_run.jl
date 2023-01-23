module ANOPP2Run

using AcousticMetrics: AcousticMetrics
import AcousticAnalogies, DelimitedFiles, FLOWMath
using StaticArrays: @SVector, SVector
using KinematicCoordinateTransformations: KinematicCoordinateTransformations, compose
using Printf: @sprintf
using LinearAlgebra: ×

include(joinpath(@__DIR__, "gen_test_data", "gen_ccblade_data", "constants.jl"))
using .CCBladeTestCaseConstants
ccbc = CCBladeTestCaseConstants

# This is a function that will do the AcousticAnalogies.jl noise calculation with the new
# observer stuff.
function cf1a_noise(num_blades, v, omega, radii, dradii, cs_area, fn, fc, stationary_observer, theta, f_interp)
    t0 = 0.0
    rot_axis = @SVector [0.0, 0.0, 1.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    x0 = SVector{3}([cos(theta), 0.0, sin(theta)].*100.0.*12.0.*0.0254)  # 100 ft in meters
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    v0_hub = SVector{3}(v.*rot_axis)
    num_src_times = 256
    num_obs_times = 2*num_src_times

    # Blade Passing Period.
    bpp = 2*pi/omega/num_blades
    src_time_range = 5.0*bpp
    obs_time_range = 4.0*bpp

    if stationary_observer
        obs = AcousticAnalogies.StationaryAcousticObserver(x0)
    else
        obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0, x0, v0_hub)
    end

    rot_trans = KinematicCoordinateTransformations.SteadyRotXTransformation(t0, omega, 0.0)
    global_trans = KinematicCoordinateTransformations.ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))
    const_vel_trans = KinematicCoordinateTransformations.ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # This is just an array of the angular offsets of each blade.
    θs = 2*pi/num_blades.*(0:(num_blades-1))

    dt = src_time_range/(num_src_times - 1)
    src_times = t0 .+ (0:num_src_times-1).*dt

    θs = reshape(θs, 1, 1, :)
    radii = reshape(radii, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    cs_area = reshape(cs_area, 1, :, 1)
    fn = reshape(fn, 1, :, 1)
    fc = reshape(fc, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This isn't really necessary.

    # Get all the transformations!
    trans = compose.(src_times, Ref(const_vel_trans), compose.(src_times, Ref(global_trans), Ref(rot_trans)))

    # Transform the source elements.
    ses = AcousticAnalogies.CompactSourceElement.(ccbc.rho, ccbc.c0, radii, θs, dradii, cs_area, -fn, 0.0, fc, src_times) .|> trans

    # Do the acoustics.
    apth = AcousticAnalogies.f1a.(ses, Ref(obs))

    # Combine all the acoustic pressure time histories into one.
    apth_total = AcousticAnalogies.combine(apth, obs_time_range, num_obs_times, 1; f_interp=f_interp)

    return AcousticMetrics.time(apth_total), AcousticAnalogies.pressure_monopole(apth_total), AcousticAnalogies.pressure_dipole(apth_total)
end

function get_results(; stationary_observer, theta, f_interp, rpm, irpm)
    apth_ylims_xrotor = [
        (-0.0001, 0.0001), (-0.0005, 0.0005), (-0.0020, 0.0020),
        (-0.005, 0.005), (-0.010, 0.010), (-0.020, 0.020), (-0.05, 0.05),
        (-0.10, 0.10), (-0.20, 0.20), (-0.5, 0.5), (-1.0, 1.0)]
    dradii = AcousticAnalogies.get_dradii(ccbc.radii, ccbc.Rhub, ccbc.Rtip)
    cs_area =  ccbc.area_over_chord_squared .* ccbc.chord.^2

    omega = rpm*(2*pi/60.0)

    # Get the normal and circumferential loading from the CCBlade output.
    fname = joinpath(@__DIR__, "gen_test_data", "gen_ccblade_data", "ccblade_omega$(@sprintf "%02d" irpm).csv")
    data = DelimitedFiles.readdlm(fname, ',')
    fn = data[:, 1]
    fc = data[:, 2]

    # Blade passing period.
    bpp = 2*pi/omega/ccbc.num_blades
    
    # Calculate the noise with AcousticAnalogies.jl.
    obs_time, p_thickness, p_loading = cf1a_noise(ccbc.num_blades, ccbc.v, omega, ccbc.radii, dradii, cs_area, fn, fc, stationary_observer, theta, f_interp)
    t0 = obs_time[1]

    # Nondimensionalize the observer time with the blade passing period.
    obs_time = (obs_time .- t0)./bpp

    # Read the ANOPP2 data.
    if theta ≈ 0.0
        if stationary_observer
            fname = joinpath(@__DIR__, "gen_test_data", "anopp2_omega$(@sprintf "%02d" irpm).csv")
            data = DelimitedFiles.readdlm(fname, ',')
        else
            fname = joinpath(@__DIR__, "gen_test_data", "anopp2_const_vel_omega$(@sprintf "%02d" irpm).csv")
            data = DelimitedFiles.readdlm(fname, ',')
        end
    else
        theta_int = Int(round(theta*180.0/pi))
        if stationary_observer
            fname = joinpath(@__DIR__, "gen_test_data", "anopp2_omega$(@sprintf "%02d" irpm)_theta$(theta_int).csv")
            data = DelimitedFiles.readdlm(fname, ',')
        else
            fname = joinpath(@__DIR__, "gen_test_data", "anopp2_const_vel_omega$(@sprintf "%02d" irpm)_theta$(theta_int).csv")
            data = DelimitedFiles.readdlm(fname, ',')
        end
    end
    t_a2 = data[:, 1]
    p_monopole_a2 = data[:, 2]
    p_dipole_a2 = data[:, 3]

    # Nondimensionalize the observer time with the blade passing period.
    t_a2 = (t_a2 .- t0)./bpp

    # Let's interpolate both the AcousticAnalogies.jl and ANOPP2 data onto a common time
    # domain.
    t = range(0.0, 1.0, length=128)
    p_thickness_interp = FLOWMath.akima(obs_time, p_thickness, t)
    p_loading_interp = FLOWMath.akima(obs_time, p_loading, t)
    p_monopole_a2_interp = FLOWMath.akima(t_a2, p_monopole_a2, t)
    p_dipole_a2_interp = FLOWMath.akima(t_a2, p_dipole_a2, t)

    return t, p_thickness_interp, p_loading_interp, p_monopole_a2_interp, p_dipole_a2_interp
end

end # module
