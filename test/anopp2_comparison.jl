module ANOPP2Comparison

# const DO_PLOTS = false

using AcousticMetrics
import AcousticAnalogies, DelimitedFiles, FLOWMath, Test
using KinematicCoordinateTransformations
using LinearAlgebra: ×
using StaticArrays
# if DO_PLOTS
#     import Plots
# end
using Printf: @sprintf

include(joinpath(@__DIR__, "gen_test_data", "gen_ccblade_data", "constants.jl"))
using .CCBladeTestCaseConstants
ccbc = CCBladeTestCaseConstants

Test.@testset "ANOPP2 Comparison" begin

    function omega_sweep_with_acoustics(; stationary_observer, theta, f_interp)
		apth_ylims_xrotor = [
			(-0.0001, 0.0001), (-0.0005, 0.0005), (-0.0020, 0.0020),
			(-0.005, 0.005), (-0.010, 0.010), (-0.020, 0.020), (-0.05, 0.05),
			(-0.10, 0.10), (-0.20, 0.20), (-0.5, 0.5), (-1.0, 1.0)]
        dradii = AcousticAnalogies.get_dradii(ccbc.radii, ccbc.Rhub, ccbc.Rtip)
        cs_area =  ccbc.area_over_chord_squared .* ccbc.chord.^2

        rpm = 200.0:200.0:2200.0  # rev/min

        for i in eachindex(rpm)
            omega = rpm[i]*(2*pi/60.0)

            # Get the normal and circumferential loading from the CCBlade output.
            data = DelimitedFiles.readdlm("gen_test_data/gen_ccblade_data/ccblade_omega$(@sprintf "%02d" i).csv", ',')
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
                    data = DelimitedFiles.readdlm("gen_test_data/anopp2_omega$(@sprintf "%02d" i).csv", ',')
                else
                    data = DelimitedFiles.readdlm("gen_test_data/anopp2_const_vel_omega$(@sprintf "%02d" i).csv", ',')
                end
            else
                theta_int = Int(round(theta*180.0/pi))
                if stationary_observer
                    data = DelimitedFiles.readdlm("gen_test_data/anopp2_omega$(@sprintf "%02d" i)_theta$(theta_int).csv", ',')
                else
                    data = DelimitedFiles.readdlm("gen_test_data/anopp2_const_vel_omega$(@sprintf "%02d" i)_theta$(theta_int).csv", ',')
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

            # Now compare the results from the two codes.
            max_aerr_thickness = maximum(abs.(p_thickness_interp - p_monopole_a2_interp))
            p_thickness_ref = maximum(p_monopole_a2_interp) - minimum(p_monopole_a2_interp)
            max_aerr_thickness_scaled = max_aerr_thickness / p_thickness_ref
            tol = 0.01
            Test.@test max_aerr_thickness_scaled < tol
            if ! (max_aerr_thickness_scaled < tol)
                println("stationary_observer = $(stationary_observer), theta = $(theta*180.0/pi) deg, rpm = $(rpm[i]), thickness aerr_scaled = $(max_aerr_thickness_scaled)")
            end

            p_loading_ref = maximum(p_dipole_a2_interp) - minimum(p_dipole_a2_interp)
            max_aerr_loading = maximum(abs.(p_loading_interp - p_dipole_a2_interp))
            max_aerr_loading_scaled = max_aerr_loading / p_loading_ref
            tol = 0.01
            if theta ≈ 0.0
                tol = 0.01
            else
                if rpm[i] ≈ 200.0
                    # The 200.0 RPM case is strange. The thrust and torque
                    # are actually negative, and the loading noise looks
                    # especially different from what ANOPP2 predicts.
                    tol = 0.045
                else
                    tol = 0.01
                end
            end
            Test.@test max_aerr_loading_scaled < tol
            if ! (max_aerr_loading_scaled < tol)
                println("stationary_observer = $(stationary_observer), theta = $(theta*180.0/pi) deg, rpm = $(rpm[i]), loading aerr_scaled = $(max_aerr_loading_scaled)")
            end

            # if DO_PLOTS
            #     # Create a plot.
            #     p_apth = Plots.plot(legend=:topleft, foreground_color_legend=nothing, background_color_legend=nothing)

            #     # Plot original data.
            #     Plots.plot!(p_apth, obs_time, p_thickness, label="thickness", seriescolor=:blue, markershape=:x)
            #     Plots.plot!(p_apth, obs_time, p_loading, label="loading", seriescolor=:red, markershape=:x)
            #     Plots.plot!(p_apth, t_a2, p_monopole_a2, label="ANOPP2 monopole", seriescolor=:blue, markershape=:+, markerstrokecolor=:blue)
            #     Plots.plot!(p_apth, t_a2, p_dipole_a2, label="ANOPP2 dipole", seriescolor=:red, markershape=:+, markerstrokecolor=:red)
            #     # Plot interpolated data.
            #     Plots.plot!(p_apth, t, p_thickness_interp, label="thickness, interp", seriescolor=:blue, markershape=:star4, markerstrokecolor=:blue, markercolor=nothing)
            #     Plots.plot!(p_apth, t, p_loading_interp, label="loading, interp", seriescolor=:red, markershape=:star4, markerstrokecolor=:red, markercolor=nothing)
            #     Plots.plot!(p_apth, t, p_monopole_a2_interp, label="ANOPP2 monopole, interp", seriescolor=:blue, markershape=:star6, markerstrokecolor=:blue, markercolor=nothing)
            #     Plots.plot!(p_apth, t, p_dipole_a2_interp, label="ANOPP2 dipole, interp", seriescolor=:red, markershape=:star6, markerstrokecolor=:red, markercolor=nothing)

            #     # Make the plot look good.
            #     if theta ≈ 0.0
            #         Plots.ylims!(p_apth, apth_ylims_xrotor[i])
            #     end
            #     Plots.xlims!(p_apth, (0.0, 1.0))
            #     Plots.xlabel!(p_apth, "t/t_blade_pass")
            #     Plots.ylabel!(p_apth, "acoustic pressure, Pa")

            #     # Save the plot.
            #     theta_int = Int(round(theta*180.0/pi))
            #     if stationary_observer
            #         if f_interp === FLOWMath.akima
            #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_akima.png")
            #         elseif f_interp === FLOWMath.linear
            #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_linear.png")
            #         else
            #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_other.png")
            #         end
            #     else
            #         if f_interp === FLOWMath.akima
            #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_akima.png")
            #         elseif f_interp === FLOWMath.linear
            #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_linear.png")
            #         else
            #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[i])))rpm_theta$(theta_int)_other.png")
            #         end
            #     end
            # end
        end 

    end

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

        rot_trans = SteadyRotXTransformation(t0, omega, 0.0)
        global_trans = ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))
        const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

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

    # Actually run the tests.
    for f_interp in [FLOWMath.akima, FLOWMath.linear]
        for theta in [0.0, -45.0*pi/180.0]
            for stationary_observer in [true, false]
                omega_sweep_with_acoustics(stationary_observer=stationary_observer, theta=theta, f_interp=f_interp)
            end
        end
    end

end

end # module
