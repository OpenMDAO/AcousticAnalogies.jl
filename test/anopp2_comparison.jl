module ANOPP2Comparison

# const DO_PLOTS = false

using FLOWMath: FLOWMath
using Test: Test
# if DO_PLOTS
#     import Plots
# end
using Printf: @sprintf

include(joinpath(@__DIR__, "anopp2_run.jl"))
using .ANOPP2Run: get_results

Test.@testset "ANOPP2 Comparison" begin

    function compare_results(; stationary_observer, theta, f_interp, rpm, irpm)
        t, p_thickness_interp, p_loading_interp, p_monopole_a2_interp, p_dipole_a2_interp = get_results(; stationary_observer, theta, f_interp, rpm, irpm)

        # Now compare the results from the two codes.
        max_aerr_thickness = maximum(abs.(p_thickness_interp - p_monopole_a2_interp))
        p_thickness_ref = maximum(p_monopole_a2_interp) - minimum(p_monopole_a2_interp)
        max_aerr_thickness_scaled = max_aerr_thickness / p_thickness_ref
        tol = 0.01
        Test.@test max_aerr_thickness_scaled < tol
        if ! (max_aerr_thickness_scaled < tol)
            println("stationary_observer = $(stationary_observer), theta = $(theta*180.0/pi) deg, rpm = $(rpm), thickness aerr_scaled = $(max_aerr_thickness_scaled)")
        end

        p_loading_ref = maximum(p_dipole_a2_interp) - minimum(p_dipole_a2_interp)
        max_aerr_loading = maximum(abs.(p_loading_interp - p_dipole_a2_interp))
        max_aerr_loading_scaled = max_aerr_loading / p_loading_ref
        tol = 0.01
        if theta ≈ 0.0
            tol = 0.01
        else
            if rpm ≈ 200.0
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
            println("stationary_observer = $(stationary_observer), theta = $(theta*180.0/pi) deg, rpm = $(rpm), loading aerr_scaled = $(max_aerr_loading_scaled)")
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
        #         Plots.ylims!(p_apth, apth_ylims_xrotor[irpm])
        #     end
        #     Plots.xlims!(p_apth, (0.0, 1.0))
        #     Plots.xlabel!(p_apth, "t/t_blade_pass")
        #     Plots.ylabel!(p_apth, "acoustic pressure, Pa")

        #     # Save the plot.
        #     theta_int = Int(round(theta*180.0/pi))
        #     if stationary_observer
        #         if f_interp === FLOWMath.akima
        #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_akima.png")
        #         elseif f_interp === FLOWMath.linear
        #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_linear.png")
        #         else
        #             Plots.savefig(p_apth, "p_apth_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_other.png")
        #         end
        #     else
        #         if f_interp === FLOWMath.akima
        #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_akima.png")
        #         elseif f_interp === FLOWMath.linear
        #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_linear.png")
        #         else
        #             Plots.savefig(p_apth, "p_apth_const_vel_$(@sprintf "%04d" Int(round(rpm[irpm])))rpm_theta$(theta_int)_other.png")
        #         end
        #     end
        # end
    end

    # Actually run the tests.
    for f_interp in [FLOWMath.akima, FLOWMath.linear]
        for theta in [0.0, -45.0*pi/180.0]
            for stationary_observer in [true, false]
                for (i, rpm) in enumerate(200.0:200.0:2200.0)  # rev/min
                    # get_results(stationary_observer=stationary_observer, theta=theta, f_interp=f_interp, rpm=rpm, irpm=i)
                    compare_results(; stationary_observer=stationary_observer, theta=theta, f_interp=f_interp, rpm=rpm, irpm=i)
                end
            end
        end
    end

end

end # module
