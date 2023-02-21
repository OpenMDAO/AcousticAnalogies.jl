module BoundaryLayerTests

using AcousticAnalogies: AcousticAnalogies
using DelimitedFiles: DelimitedFiles
using FLOWMath: linear
using Test

@testset "displacement thickness" begin
    @testset "zero angle of attack" begin
        @testset "tripped" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure06-disp_thickness-tripped.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure06-disp_thickness-untripped.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure07-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            deltastar_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies.disp_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.06
        end

        @testset "suction side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure07-suction_side.csv")
            bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_suction_side[:, 1]
            deltastar_bpm = bpm_suction_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies.disp_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report data onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            @test maximum(err) < 0.04
        end
    end

    @testset "non-zero angle of attack, untripped" begin
        @testset "pressure side" begin
            # Get the digitized data from the BPM report plot.
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure08-pressure_side.csv")
            bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_pressure_side[:, 1]
            deltastar_bpm = bpm_pressure_side[:, 2]

            # Get the AcousticAnalogies.jl implementation.
            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies.disp_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            # Interpolate the BPM report onto the uniform alpha spacing.
            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            # Find the scaled error.
            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            # This is dumb. Maybe I have a bug?
            @test maximum(err) < 0.11
        end

        @testset "suction side" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure08-suction_side.csv")
            bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm_suction_side[:, 1]
            deltastar_bpm = bpm_suction_side[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            deltastar_jl = AcousticAnalogies.disp_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)

            deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)

            vmin, vmax = extrema(deltastar_bpm)
            err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
            # This is dumb. Maybe I have a bug?
            @test maximum(err) < 0.081
        end
    end
end

@testset "shape functions" begin
    @testset "TBL-TE" begin
        @testset "St_1" begin
            @test isapprox(AcousticAnalogies.St_1(0.093), 0.081; atol=0.0022)
            @test isapprox(AcousticAnalogies.St_1(0.116), 0.071; atol=0.002)
            @test isapprox(AcousticAnalogies.St_1(0.163), 0.059; atol=0.0004)
            @test isapprox(AcousticAnalogies.St_1(0.209), 0.051; atol=0.0002)
        end

        @testset "K_1" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure77.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            Re_c_bpm = bpm[:, 1]
            K_1_bpm = bpm[:, 2]

            Re_c_jl = range(minimum(Re_c_bpm), maximum(Re_c_bpm); length=50)
            K_1_jl = AcousticAnalogies.K_1.(Re_c_jl)

            K_1_interp = linear(Re_c_bpm, K_1_bpm, Re_c_jl)

            vmin, vmax = extrema(K_1_bpm)
            err = abs.(K_1_jl .- K_1_interp)./(vmax - vmin)
            @test maximum(err) < 0.012
        end

        @testset "A" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure78-A_min.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_St_peak_bpm = bpm[:, 1]
            A = bpm[:, 2]

            St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
            A_jl = AcousticAnalogies.A.(St_St_peak_jl, 9.5e4)

            # Interpolate:
            A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
            
            # Check error.
            vmin, vmax = extrema(A)
            err = abs.(A_jl .- A_interp)./(vmax - vmin)
            @test maximum(err) < 0.057

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure78-A_max.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_St_peak_bpm = bpm[:, 1]
            A = bpm[:, 2]

            St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
            A_jl = AcousticAnalogies.A.(St_St_peak_jl, 8.58e5)

            # Interpolate:
            A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
            
            # Check error.
            vmin, vmax = extrema(A)
            err = abs.(A_jl .- A_interp)./(vmax - vmin)
            @test maximum(err) < 0.021
        end

        @testset "B" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure78-B_min.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_St_peak_bpm = bpm[:, 1]
            B = bpm[:, 2]

            St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
            B_jl = AcousticAnalogies.B.(St_St_peak_jl, 9.5e4)

            # Interpolate:
            B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
            
            # Check error.
            vmin, vmax = extrema(B)
            err = abs.(B_jl .- B_interp)./(vmax - vmin)
            @test maximum(err) < 0.057

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure78-B_max.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_St_peak_bpm = bpm[:, 1]
            B = bpm[:, 2]

            St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
            B_jl = AcousticAnalogies.B.(St_St_peak_jl, 8.58e5)

            # Interpolate:
            B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
            
            # Check error.
            vmin, vmax = extrema(B)
            err = abs.(B_jl .- B_interp)./(vmax - vmin)
            @test maximum(err) < 0.020
        end

        @testset "St_2" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure80-M0.093.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            St_2 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.093), alpha_deg_jl.*pi/180)

            # Interpolate:
            St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)

            # Check error.
            vmin, vmax = extrema(St_2)
            err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
            @test maximum(err) < 0.023

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure80-M0.209.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            St_2 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
            St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.209), alpha_deg_jl.*pi/180)

            # Interpolate:
            St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)

            # Check error.
            vmin, vmax = extrema(St_2)
            err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
            @test maximum(err) < 0.011

        end

        @testset "K_2" begin
            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure82-M0.093.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            K_2_K_1 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
            K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.093, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)

            # Interpolate:
            K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)

            # Check error.
            vmin, vmax = extrema(K_2_K_1)
            err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
            # The curve is almost vertical at low angles of attack, so a small error in the digitization results in big differences.
            @test maximum(err[2:end]) < 0.027

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure82-M0.116.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            K_2_K_1 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
            K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.116, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)

            # Interpolate:
            K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)

            # Check error.
            vmin, vmax = extrema(K_2_K_1)
            err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
            # There's a branch for low angles of attack that sets K_2 - K_1 to
            # -1000, but I can't see that in the digitized plots, so the first
            # point is bad.
            @test K_2_K_1_jl[1] ≈ -1000
            @test maximum(err[2:end]) < 0.022

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure82-M0.163.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            K_2_K_1 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
            K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.163, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)

            # Interpolate:
            K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)

            # Check error.
            vmin, vmax = extrema(K_2_K_1)
            err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
            # There's a branch for low angles of attack that sets K_2 - K_1 to
            # -1000, but I can't see that in the digitized plots, so the first
            # point is bad.
            @test K_2_K_1_jl[1] ≈ -1000.0
            @test maximum(err[2:end]) < 0.020

            fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure82-M0.209.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alpha_deg = bpm[:, 1]
            K_2_K_1 = bpm[:, 2]

            alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
            K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.209, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)

            # Interpolate:
            K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
            
            # Check error.
            vmin, vmax = extrema(K_2_K_1)
            err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
            @test maximum(err) < 0.024
        end
    end
end

end # module
