module BPMShapeFunctionTests

using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: DelimitedFiles
using FLOWMath: linear
using Test

@testset "shape functions" begin

    @testset "TBL-TE" begin
        @testset "St_1" begin
            @test isapprox(AcousticAnalogies.St_1(0.093), 0.081; atol=0.0022)
            @test isapprox(AcousticAnalogies.St_1(0.116), 0.071; atol=0.002)
            @test isapprox(AcousticAnalogies.St_1(0.163), 0.059; atol=0.0004)
            @test isapprox(AcousticAnalogies.St_1(0.209), 0.051; atol=0.0002)
        end

        @testset "K_1" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure77.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-A_min.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-A_max.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-B_min.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-B_max.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure80-M0.093.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure80-M0.209.csv")
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
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.093.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.116.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.163.csv")
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

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.209.csv")
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

    @testset "LBL-VS" begin
        @testset "St_1_prime" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure86-St_1_prime.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            Re_c_bpm = bpm[:, 1]
            St_1_prime_bpm = bpm[:, 2]

            Re_c_jl = 10.0.^(range(4, 7; length=100))
            St_1_prime_jl = AcousticAnalogies.St_1_prime.(Re_c_jl)

            St_1_prime_interp = linear(Re_c_bpm, St_1_prime_bpm, Re_c_jl)
            vmin, vmax = extrema(St_1_prime_bpm)
            err = abs.(St_1_prime_interp .- St_1_prime_jl)./(vmax - vmin)
            @test maximum(err) < 0.04
        end

        @testset "G1" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure85-G1.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            e_bpm = bpm[:, 1]
            G1_bpm = bpm[:, 2]

            e_jl = 10.0.^(range(-1, 1; length=101))
            G1_jl = AcousticAnalogies.G1.(e_jl)

            G1_interp = linear(e_jl, G1_jl, e_bpm)
            vmin, vmax = extrema(G1_bpm)
            err = abs.(G1_interp .- G1_bpm)./(vmax - vmin)
            @test maximum(err) < 0.033
        end

        @testset "St_peak_prime_alphastar" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure87.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            alphastar_bpm = bpm[:, 1]
            St_peak_ratio_bpm = bpm[:, 2]

            St_1_prime = 0.25  # Just make up a value, since we're multiplying and then dividing by it anyway.
            alphastar_jl = range(0.0*pi/180, 7.0*pi/180; length=21)
            St_peak_ratio_jl = AcousticAnalogies.St_peak_prime.(St_1_prime, alphastar_jl)./St_1_prime

            St_peak_ratio_interp = linear(alphastar_jl.*180/pi, St_peak_ratio_jl, alphastar_bpm)
            vmin, vmax = extrema(St_peak_ratio_bpm)
            err = abs.(St_peak_ratio_interp .- St_peak_ratio_bpm)./(vmax - vmin)
            @test maximum(err) < 0.031
        end

        @testset "G2" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure89.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            Re_ratio_bpm = bpm[:, 1]
            G2_bpm = bpm[:, 2]

            Re_ratio_jl = 10.0.^range(-1, 1, length=51)
            G2_jl = AcousticAnalogies.G2.(Re_ratio_jl)

            G2_interp = linear(Re_ratio_jl, G2_jl, Re_ratio_bpm)
            vmin, vmax = extrema(G2_interp)
            err = abs.(G2_interp .- G2_bpm)./(vmax - vmin)
            @test maximum(err) < 0.024
        end

        @testset "G2 + G3" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure88-G2-alpha0.csv")
            alphastar = 0.0*pi/180
            bpm = DelimitedFiles.readdlm(fname, ',')
            Re_c_bpm = bpm[:, 1]
            G2_bpm = bpm[:, 2]

            Re_c_jl = 10.0.^range(log10(first(Re_c_bpm)), log10(last(Re_c_bpm)), length=51)
            Re_c0 = AcousticAnalogies.Re_c0(alphastar)
            Re_ratio_jl = Re_c_jl./Re_c0
            G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ AcousticAnalogies.G3.(alphastar)

            G2_interp = linear(Re_c_jl, G2_jl, Re_c_bpm)
            vmin, vmax = extrema(G2_interp)
            err = abs.(G2_interp .- G2_bpm)./(vmax - vmin)
            @test maximum(err) < 0.013

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure88-G2-alpha6.csv")
            alphastar = 6.0*pi/180
            bpm = DelimitedFiles.readdlm(fname, ',')
            Re_c_bpm = bpm[:, 1]
            G2_bpm = bpm[:, 2]

            Re_c_jl = 10.0.^range(log10(first(Re_c_bpm)), log10(last(Re_c_bpm)), length=51)
            Re_c0 = AcousticAnalogies.Re_c0(alphastar)
            Re_ratio_jl = Re_c_jl./Re_c0
            G2_jl = AcousticAnalogies.G2.(Re_ratio_jl) .+ AcousticAnalogies.G3.(alphastar)

            G2_interp = linear(Re_c_jl, G2_jl, Re_c_bpm)
            vmin, vmax = extrema(G2_interp)
            err = abs.(G2_interp .- G2_bpm)./(vmax - vmin)
            @test maximum(err) < 0.030
        end

    end

    @testset "TEB-VS" begin
        @testset "St_3prime_peak" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure95-0Psi.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            h_over_deltastar_0Psi = bpm[:, 1]
            St_3prime_peak_0Psi = bpm[:, 2]
             
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure95-14Psi.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            h_over_deltastar_14Psi = bpm[:, 1]
            St_3prime_peak_14Psi = bpm[:, 2]

            h_over_deltastar_jl = 10.0.^(range(-1, 1; length=101))
            St_3prime_peak_0Psi_jl = AcousticAnalogies.St_3prime_peak.(h_over_deltastar_jl, 0.0*pi/180)
            St_3prime_peak_14Psi_jl = AcousticAnalogies.St_3prime_peak.(h_over_deltastar_jl, 14.0*pi/180)

            St_3prime_peak_0Psi_interp = linear(h_over_deltastar_jl, St_3prime_peak_0Psi_jl, h_over_deltastar_0Psi)
            vmin, vmax = extrema(St_3prime_peak_0Psi)
            err = abs.(St_3prime_peak_0Psi_interp .- St_3prime_peak_0Psi)./(vmax - vmin)
            @test maximum(err) < 0.070

            St_3prime_peak_14Psi_interp = linear(h_over_deltastar_jl, St_3prime_peak_14Psi_jl, h_over_deltastar_14Psi)
            vmin, vmax = extrema(St_3prime_peak_14Psi)
            err = abs.(St_3prime_peak_14Psi_interp .- St_3prime_peak_14Psi)./(vmax - vmin)
            @test maximum(err) < 0.049
        end

        @testset "G4" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure96-0Psi.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            h_over_deltastar_0Psi = bpm[:, 1]
            G4_0Psi = bpm[:, 2]
             
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure96-14Psi.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            h_over_deltastar_14Psi = bpm[:, 1]
            G4_14Psi = bpm[:, 2]

            h_over_deltastar_jl = 10.0.^(range(-1, 1; length=51))
            G4_0Psi_jl = AcousticAnalogies.G4.(h_over_deltastar_jl, 0.0*pi/180)
            G4_14Psi_jl = AcousticAnalogies.G4.(h_over_deltastar_jl, 14.0*pi/180)

            G4_0Psi_interp = linear(h_over_deltastar_jl, G4_0Psi_jl, h_over_deltastar_0Psi)
            vmin, vmax = extrema(G4_0Psi)
            err = abs.(G4_0Psi_interp .- G4_0Psi)./(vmax - vmin)
            @test maximum(err) < 0.033

            G4_14Psi_interp = linear(h_over_deltastar_jl, G4_14Psi_jl, h_over_deltastar_14Psi)
            vmin, vmax = extrema(G4_14Psi)
            err = abs.(G4_14Psi_interp .- G4_14Psi)./(vmax - vmin)
            @test maximum(err) < 0.024

        end

        @testset "G5, Psi = 14°" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar0p25.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p25 = bpm[:, 1]
            G5_14Psi_h_over_deltastar_avg0p25 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar0p43.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p43 = bpm[:, 1]
            G5_14Psi_h_over_deltastar_avg0p43 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar0p50.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p50 = bpm[:, 1]
            G5_14Psi_h_over_deltastar_avg0p50 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar0p54.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p54 = bpm[:, 1]
            G5_14Psi_h_over_deltastar_avg0p54 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar0p62.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p62 = bpm[:, 1]
            G5_14Psi_h_over_deltastar_avg0p62 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi14-h_over_deltastar1p20.csv")
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

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p25_jl, St_3prime_over_St_3prime_peak_0p25)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg0p25)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg0p25)/(vmax - vmin)
            @test maximum(err) < 0.074

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p43_jl, St_3prime_over_St_3prime_peak_0p43)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg0p43)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg0p43)/(vmax - vmin)
            @test maximum(err) < 0.072

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p50_jl, St_3prime_over_St_3prime_peak_0p50)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg0p50)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg0p50)/(vmax - vmin)
            @test maximum(err) < 0.072

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p54_jl, St_3prime_over_St_3prime_peak_0p54)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg0p54)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg0p54)/(vmax - vmin)
            @test maximum(err) < 0.074

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg0p62_jl, St_3prime_over_St_3prime_peak_0p62)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg0p62)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg0p62)/(vmax - vmin)
            @test maximum(err) < 0.073

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_14Psi_h_over_deltastar_avg1p20_jl, St_3prime_over_St_3prime_peak_1p20)
            vmin, vmax = extrema(G5_14Psi_h_over_deltastar_avg1p20)
            err = abs.(interp .- G5_14Psi_h_over_deltastar_avg1p20)/(vmax - vmin)
            # The lower end of this case is really bad.
            # Not sure why. :-(
            @test maximum(err[1:22]) < 0.31
            @test maximum(err[23:end]) < 0.087
        end

        @testset "G5, Psi = 0°" begin
            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar0p25.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p25 = bpm[:, 1]
            G5_0Psi_h_over_deltastar_avg0p25 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar0p43.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p43 = bpm[:, 1]
            G5_0Psi_h_over_deltastar_avg0p43 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar0p50.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p50 = bpm[:, 1]
            G5_0Psi_h_over_deltastar_avg0p50 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar0p54.csv")
            bpm = DelimitedFiles.readdlm(fname, ',')
            St_3prime_over_St_3prime_peak_0p54 = bpm[:, 1]
            G5_0Psi_h_over_deltastar_avg0p54 = bpm[:, 2]

            # fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar0p62.csv")
            # bpm = DelimitedFiles.readdlm(fname, ',')
            # St_3prime_over_St_3prime_peak_0p62 = bpm[:, 1]
            # G5_0Psi_h_over_deltastar_avg0p62 = bpm[:, 2]

            fname = joinpath(@__DIR__, "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure97-Psi0-h_over_deltastar1p20.csv")
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

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p25_jl, St_3prime_over_St_3prime_peak_0p25)
            vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg0p25)
            err = abs.(interp .- G5_0Psi_h_over_deltastar_avg0p25)/(vmax - vmin)
            @test maximum(err) < 0.030

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p43_jl, St_3prime_over_St_3prime_peak_0p43)
            vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg0p43)
            err = abs.(interp .- G5_0Psi_h_over_deltastar_avg0p43)/(vmax - vmin)
            @test maximum(err) < 0.026

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p50_jl, St_3prime_over_St_3prime_peak_0p50)
            vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg0p50)
            err = abs.(interp .- G5_0Psi_h_over_deltastar_avg0p50)/(vmax - vmin)
            @test maximum(err) < 0.037

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p54_jl, St_3prime_over_St_3prime_peak_0p54)
            vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg0p54)
            err = abs.(interp .- G5_0Psi_h_over_deltastar_avg0p54)/(vmax - vmin)
            @test maximum(err) < 0.037

            # interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg0p62_jl, St_3prime_over_St_3prime_peak_0p62)
            # vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg0p62)
            # err = abs.(interp .- G5_0Psi_h_over_deltastar_avg0p62)/(vmax - vmin)
            # @test maximum(err) < 0.073

            interp = linear(St_3prime_over_St_3prime_peak_jl, G5_0Psi_h_over_deltastar_avg1p20_jl, St_3prime_over_St_3prime_peak_1p20)
            vmin, vmax = extrema(G5_0Psi_h_over_deltastar_avg1p20)
            err = abs.(interp .- G5_0Psi_h_over_deltastar_avg1p20)/(vmax - vmin)
            @test maximum(err) < 0.045
        end

    end

end

end # module
