module BoundaryLayerTests

using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: DelimitedFiles
using FLOWMath: linear
using Test

@testset "boundary layer thickness" begin
    @testset "zero angle of attack" begin
        @testset "untripped" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-bl_thickness-untripped.csv")
            bpm_untripped = DelimitedFiles.readdlm(fname, ',')
            Re_c_1e6 = bpm_untripped[:, 1]
            delta0_c = bpm_untripped[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
            delta0_c_jl = AcousticAnalogies.bl_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
            # Interpolate the BPM report data onto the uniform Re spacing.
            delta0_c_interp = linear(Re_c_1e6, delta0_c, Re_c_1e6_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta0_c)
            err = abs.(delta0_c_jl .- delta0_c_interp)/(vmax - vmin)
            @test maximum(err) < 0.041

        end

        @testset "tripped" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-bl_thickness-tripped.csv")
            bpm_tripped = DelimitedFiles.readdlm(fname, ',')
            Re_c_1e6 = bpm_tripped[:, 1]
            delta0_c = bpm_tripped[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
            delta0_c_jl = AcousticAnalogies.bl_thickness_0.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
            # Interpolate the BPM report data onto the uniform Re spacing.
            delta0_c_interp = linear(Re_c_1e6, delta0_c, Re_c_1e6_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta0_c)
            err = abs.(delta0_c_jl .- delta0_c_interp)/(vmax - vmin)
            @test maximum(err) < 0.0086

        end
    end

    @testset "non-zero angle of attack, tripped" begin
        @testset "pressure side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-bl_thickness-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            delta_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            delta_jl = AcousticAnalogies._bl_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta_bpm)
            err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.032
        end

        @testset "suction side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-bl_thickness-suction_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            delta_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            delta_jl = AcousticAnalogies._bl_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta_bpm)
            err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.023
        end
    end

    @testset "non-zero angle of attack, untripped" begin
        @testset "pressure side" begin
            # Non-zero boundary layer thickness for the pressure side is the same for tripped and untripped, so getting the data for the tripped case.
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-bl_thickness-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            delta_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            delta_jl = AcousticAnalogies._bl_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta_bpm)
            err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.032
        end

        @testset "suction side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-bl_thickness-suction_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            delta_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            delta_jl = AcousticAnalogies._bl_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta_bpm)
            err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.028
        end
    end

    @testset "positive/negative angle of attack" begin
        for bl in [AcousticAnalogies.TrippedN0012BoundaryLayer(), AcousticAnalogies.UntrippedN0012BoundaryLayer()]
            for Re_c in (range(0.04, 3.0; length=30)) .* 10^6
                for alphastar_deg in 0:30
                    alphastar = alphastar_deg*pi/180
                    # For a positive angle of attack, the pressure side should be the bottom side.
                    delta_p = AcousticAnalogies.bl_thickness_p(bl, Re_c, alphastar)
                    delta_bot = AcousticAnalogies.bl_thickness_bot(bl, Re_c, alphastar)
                    @test delta_p ≈ delta_bot

                    # For a positive angle of attack, the suction side should be the top side.
                    delta_s = AcousticAnalogies.bl_thickness_s(bl, Re_c, alphastar)
                    delta_top = AcousticAnalogies.bl_thickness_top(bl, Re_c, alphastar)
                    @test delta_s ≈ delta_top

                    # But if we switch the sign on alpha, the top and bottom switch too.
                    delta_bot_neg = AcousticAnalogies.bl_thickness_bot(bl, Re_c, -alphastar)
                    @test delta_bot_neg ≈ delta_s

                    delta_top_neg = AcousticAnalogies.bl_thickness_top(bl, Re_c, -alphastar)
                    @test delta_top_neg ≈ delta_p

                    # But the value of the pressure and suction sides should never change.
                    delta_p_neg = AcousticAnalogies.bl_thickness_p(bl, Re_c, -alphastar)
                    @test delta_p_neg ≈ delta_p
                    delta_s_neg = AcousticAnalogies.bl_thickness_s(bl, Re_c, -alphastar)
                    @test delta_s_neg ≈ delta_s
                end
            end
        end
    end

    @testset "ITrip1N0012BoundaryLayer" begin

        bl = AcousticAnalogies.ITrip1N0012BoundaryLayer()
        bl_untripped = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        bl_tripped = AcousticAnalogies.TrippedN0012BoundaryLayer()
        for Re_c in (range(0.04, 3.0; length=30)) .* 10^6
            for alphastar in (-30:30) .* (pi/180)
                # Should use the untripped pressure-side and suction-side boundary layer thickness.
                @test AcousticAnalogies.bl_thickness_p(bl, Re_c, alphastar) ≈ AcousticAnalogies.bl_thickness_p(bl_untripped, Re_c, alphastar)
                @test AcousticAnalogies.bl_thickness_s(bl, Re_c, alphastar) ≈ AcousticAnalogies.bl_thickness_s(bl_untripped, Re_c, alphastar)

                # Should use the tripped displacement thickness.
                @test AcousticAnalogies.disp_thickness_p(bl, Re_c, alphastar) ≈ AcousticAnalogies.disp_thickness_p(bl_tripped, Re_c, alphastar)
                @test AcousticAnalogies.disp_thickness_s(bl, Re_c, alphastar) ≈ AcousticAnalogies.disp_thickness_s(bl_tripped, Re_c, alphastar)

            end
        end
    end

    @testset "ITrip2N0012BoundaryLayer" begin
        bl = AcousticAnalogies.ITrip2N0012BoundaryLayer()
        bl_untripped = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        bl_tripped = AcousticAnalogies.TrippedN0012BoundaryLayer()
        for Re_c in (range(0.04, 3.0; length=30)) .* 10^6
            for alphastar in (-30:30) .* (pi/180)
                # boundary layer thickness should be untripped multiplied by 0.6.
                @test AcousticAnalogies.bl_thickness_p(bl, Re_c, alphastar) ≈ 0.6*AcousticAnalogies.bl_thickness_p(bl_untripped, Re_c, alphastar)
                @test AcousticAnalogies.bl_thickness_s(bl, Re_c, alphastar) ≈ 0.6*AcousticAnalogies.bl_thickness_s(bl_untripped, Re_c, alphastar)

                # The pressure-side displacement thickness should be tripped multiplied by 0.6.
                @test AcousticAnalogies.disp_thickness_p(bl, Re_c, alphastar) ≈ 0.6*AcousticAnalogies.disp_thickness_p(bl_tripped, Re_c, alphastar)
                # The suction-side displacement thickness should be the tripped zero-alpha displacement thickness multipled by 0.6 multiplied by the untripped suction-side to zero-alpha displacement thickness ratio.
                @test AcousticAnalogies.disp_thickness_s(bl, Re_c, alphastar) ≈ AcousticAnalogies.disp_thickness_s(bl_tripped, Re_c, 0) * 0.6 * AcousticAnalogies.disp_thickness_s(bl_untripped, Re_c, alphastar) / AcousticAnalogies.disp_thickness_s(bl_untripped, Re_c, 0)
            end
        end
    end

    @testset "ITrip3N0012BoundaryLayer" begin
        bl = AcousticAnalogies.ITrip3N0012BoundaryLayer()
        bl_untripped = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        bl_tripped = AcousticAnalogies.TrippedN0012BoundaryLayer()
        for Re_c in (range(0.04, 3.0; length=30)) .* 10^6
            for alphastar in (-30:30) .* (pi/180)
                # boundary layer thickness should be untripped.
                @test AcousticAnalogies.bl_thickness_p(bl, Re_c, alphastar) ≈ AcousticAnalogies.bl_thickness_p(bl_untripped, Re_c, alphastar)
                @test AcousticAnalogies.bl_thickness_s(bl, Re_c, alphastar) ≈ AcousticAnalogies.bl_thickness_s(bl_untripped, Re_c, alphastar)

                # The pressure-side displacement thickness should be untripped multiplied by 1.48.
                @test AcousticAnalogies.disp_thickness_p(bl, Re_c, alphastar) ≈ 1.48*AcousticAnalogies.disp_thickness_p(bl_untripped, Re_c, alphastar)
                # The suction-side displacement thickness should be the untripped.
                @test AcousticAnalogies.disp_thickness_s(bl, Re_c, alphastar) ≈ AcousticAnalogies.disp_thickness_s(bl_untripped, Re_c, alphastar)
            end
        end
    end
end

@testset "displacement thickness" begin
    @testset "zero angle of attack" begin
        @testset "tripped" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-disp_thickness-tripped.csv")
            bpm_tripped = DelimitedFiles.readdlm(fname, ',')
            Re_c_1e6 = bpm_tripped[:, 1]
            deltastar0_c = bpm_tripped[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
            deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)

            # Interpolate the BPM report data onto the uniform Re spacing.
            deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar0_c)
            err = abs.(deltastar0_c_jl .- deltastar0_c_interp)/(vmax - vmin)
            @test maximum(err) < 0.05
        end

        @testset "untripped" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-disp_thickness-untripped.csv")
            bpm_untripped = DelimitedFiles.readdlm(fname, ',')
            Re_c_1e6 = bpm_untripped[:, 1]
            deltastar0_c = bpm_untripped[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
            deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
            # Interpolate the BPM report data onto the uniform Re spacing.
            deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar0_c)
            err = abs.(deltastar0_c_jl .- deltastar0_c_interp)/(vmax - vmin)
            @test maximum(err) < 0.02
        end
    end

    @testset "non-zero angle of attack, tripped" begin
        @testset "pressure side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            deltastar_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies._disp_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.06
        end

        @testset "suction side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-suction_side.csv")
            bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_suction_side[:, 1]
            deltastar_bpm = bpm_suction_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies._disp_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.04
        end
    end

    @testset "non-zero angle of attack, untripped" begin
        @testset "boundary layer thickness, pressure side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-bl_thickness-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            delta_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            delta_jl = AcousticAnalogies._bl_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report onto the uniform alpha spacing.
            delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(delta_bpm)
            err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.037
        end

        @testset "displacement thickness, pressure side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            deltastar_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies._disp_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            # This is dumb. Maybe I have a bug?
            @test maximum(err) < 0.11
        end

        @testset "displacement thinckness, suction side" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-suction_side.csv")
            bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_suction_side[:, 1]
            deltastar_bpm = bpm_suction_side[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies._disp_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            # This is dumb. Maybe I have a bug?
            @test maximum(err) < 0.081
        end
    end

    @testset "positive/negative angle of attack" begin
        for bl in [AcousticAnalogies.TrippedN0012BoundaryLayer(), AcousticAnalogies.UntrippedN0012BoundaryLayer()]
            for Re_c in (range(0.04, 3.0; length=30)) .* 10^6
                for alphastar_deg in 0:30
                    alphastar = alphastar_deg*pi/180
                    # For a positive angle of attack, the pressure side should be the bottom side.
                    deltastar_p = AcousticAnalogies.disp_thickness_p(bl, Re_c, alphastar)
                    deltastar_bot = AcousticAnalogies.disp_thickness_bot(bl, Re_c, alphastar)
                    @test deltastar_p ≈ deltastar_bot

                    # For a positive angle of attack, the suction side should be the top side.
                    deltastar_s = AcousticAnalogies.disp_thickness_s(bl, Re_c, alphastar)
                    deltastar_top = AcousticAnalogies.disp_thickness_top(bl, Re_c, alphastar)
                    @test deltastar_s ≈ deltastar_top

                    # But if we switch the sign on alpha, the top and bottom switch too.
                    deltastar_bot_neg = AcousticAnalogies.disp_thickness_bot(bl, Re_c, -alphastar)
                    @test deltastar_bot_neg ≈ deltastar_s

                    deltastar_top_neg = AcousticAnalogies.disp_thickness_top(bl, Re_c, -alphastar)
                    @test deltastar_top_neg ≈ deltastar_p

                    # But the value of the pressure and suction sides should never change.
                    deltastar_p_neg = AcousticAnalogies.disp_thickness_p(bl, Re_c, -alphastar)
                    @test deltastar_p_neg ≈ deltastar_p
                    deltastar_s_neg = AcousticAnalogies.disp_thickness_s(bl, Re_c, -alphastar)
                    @test deltastar_s_neg ≈ deltastar_s
                end
            end
        end
    end
end

end # module
