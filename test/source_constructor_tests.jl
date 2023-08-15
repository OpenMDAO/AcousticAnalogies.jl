module SourceConstructorTests

using AcousticAnalogies
using CCBlade
using DelimitedFiles
using KinematicCoordinateTransformations
using StaticArrays
using JLD2
using Test

@testset "Constructor rotation tests" begin
    @testset "CompactSourceElement" begin
        ρ0 = 1.1
        c0 = 1.2
        r = 2.0
        Δr = 0.1
        Λ = 0.2
        fn = 2.0
        fr = 3.0
        fc = 4.0
        τ = 0.1
        se_0theta = CompactSourceElement(ρ0, c0, r, 0.0, Δr, Λ, fn, fr, fc, τ)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            # Create a transformation that will undo the θ rotation.
            trans = KinematicCoordinateTransformations.SteadyRotXTransformation(τ, 0.0, -θ)
            # Create a source element with the theta rotation, then undo it.
            se = CompactSourceElement(ρ0, c0, r, θ, Δr, Λ, fn, fr, fc, τ) |> trans
            # Check that we got the same thing:
            for field in fieldnames(CompactSourceElement)
                @test getproperty(se, field) ≈ getproperty(se_0theta, field)
            end
        end
    end

    @testset "CompactSourceElement, CCBlade" begin

        # Create the CCBlade objects.
        area_per_chord2 = 0.1
        τ = 0.1
        ccblade_fname = joinpath(@__DIR__, "gen_test_data", "gen_ccblade_data", "ccblade_omega11.jld2")
        out, section, Δr, op, rotor0precone = nothing, nothing, nothing, nothing, nothing
        jldopen(ccblade_fname, "r") do f
            out = f["outs"][1]
            section = f["sections"][1]
            Δr = f["sections"][2].r - f["sections"][1].r
            op = f["ops"][1]
            rotor0precone = f["rotor"]
        end
        @test rotor0precone.precone ≈ 0.0
        Rhub = rotor0precone.Rhub
        Rtip = rotor0precone.Rtip
        num_blades = rotor0precone.B
        turbine = rotor0precone.turbine
        for positive_x_rotation in [true, false]
            se_0theta0precone = CompactSourceElement(rotor0precone, section, op, out, 0.0, Δr, area_per_chord2, τ, positive_x_rotation)
            for precone in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                rotor = CCBlade.Rotor(Rhub, Rtip, num_blades; turbine=turbine, precone=precone)
                # This is tricky: in my "normal" coordinate system, the blade is rotating around the x axis, moving axially in the positive x direction, and is initially aligned with the y axis.
                # That means that the precone should be a rotation around the negative z axis.
                # And so to undo it, we want a positive rotation around the positive z axis.
                trans_precone = SteadyRotZTransformation(τ, 0.0, precone)
                for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)
                    # Create a transformation that reverses the theta and precone rotations.
                    # The precone happens first, then theta.
                    # So to reverse it we need to do theta, then precone.
                    trans = KinematicCoordinateTransformations.compose(τ, trans_precone, trans_theta)
                    # Create a source element with the theta and precone rotations, then undo it.
                    se = CompactSourceElement(rotor, section, op, out, θ, Δr, area_per_chord2, τ, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(CompactSourceElement)
                        @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                    end
                end
            end
        end

    end

    @testset "TBLTESourceElement" begin
        # So, the way this should work: first do the twist, then do the theta rotation.
        # The twist could be either about the positive y axis or negative y axis.
        # Then the theta rotation is always about the x axis.
        c0 = 1.1
        nu = 1.2
        r = 2.0
        Δr = 0.1
        chord = 1.3
        vn = 2.0
        vr = 3.0
        vc = 4.0
        τ = 0.1
        Δτ = 0.02
        bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
        for twist_about_positive_y in [true, false]
            se_0twist0theta = TBLTESourceElement(c0, nu, r, 0.0, Δr, chord, 0.0, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
            for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

                for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    # The angle of attack depends on the twist and the fluid velocity
                    if twist_about_positive_y
                        alpha_check = ϕ - atan(-vn, -vc)
                    else
                        alpha_check = ϕ - atan(-vn, vc)
                    end
                    se = TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                    # Adjust the angles of attack to always be between -pi and pi.
                    alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                    alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                    @test alpha ≈ alpha_check

                    for field in fieldnames(TBLTESourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if field != :chord_uvec
                            @test getproperty(se, field) ≈ getproperty(se_0twist0theta, field)
                        end
                    end

                    if twist_about_positive_y
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
                        trans_phi = SteadyRotYTransformation(τ, 0.0, -ϕ)
                        chord_uvec_check = @SVector [0.0, 0.0, -1.0]
                    else
                        # If we're applying the twist about the negative y axis, then we need to do a positive rotation about the y axis to undo it.
                        trans_phi = SteadyRotYTransformation(τ, 0.0, ϕ)
                        chord_uvec_check = @SVector [0.0, 0.0, 1.0]
                    end
                    se_no_twist = se |> trans_phi
                    @test se_no_twist.chord_uvec ≈ chord_uvec_check
                end
            end
        end
    end

    @testset "TBLTESourceElement, CCBlade" begin
        # Create the CCBlade objects.
        τ = 0.1
        Δτ = 0.02
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        ccblade_fname = joinpath(@__DIR__, "gen_test_data", "gen_ccblade_data", "ccblade_omega11.jld2")
        out, section_loaded, Δr, op, rotor0precone = nothing, nothing, nothing, nothing, nothing
        jldopen(ccblade_fname, "r") do f
            out = f["outs"][1]
            section_loaded = f["sections"][1]
            Δr = f["sections"][2].r - f["sections"][1].r
            op = f["ops"][1]
            rotor0precone = f["rotor"]
            @test rotor0precone.precone ≈ 0.0
        end
        for positive_x_rotation in [true, false]
            for twist in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                section = CCBlade.Section(section_loaded.r, section_loaded.chord, twist, section_loaded.af)
                se_0theta0precone = TBLTESourceElement(rotor0precone, section, op, out, 0.0, Δr, τ, Δτ, bl, positive_x_rotation)
                for precone in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    rotor = CCBlade.Rotor(rotor0precone.Rhub, rotor0precone.Rtip, rotor0precone.B; turbine=rotor0precone.turbine, precone=precone)
                    # This is tricky: in my "normal" coordinate system, the blade is rotating around the x axis, moving axially in the positive x direction, and is initially aligned with the y axis.
                    # That means that the precone should be a rotation around the negative z axis.
                    # And so to undo it, we want a positive rotation around the positive z axis.
                    trans_precone = SteadyRotZTransformation(τ, 0.0, precone)
                    for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                        trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)
                        # Create a transformation that reverses the theta and precone rotations.
                        # The precone happens first, then theta.
                        # So to reverse it we need to do theta, then precone.
                        trans = KinematicCoordinateTransformations.compose(τ, trans_precone, trans_theta)
                        # Create a source element with the theta and precone rotations, then undo it.
                        se = TBLTESourceElement(rotor, section, op, out, θ, Δr, τ, Δτ, bl, positive_x_rotation) |> trans
                        # Check that we got the same thing:
                        for field in fieldnames(TBLTESourceElement)
                            # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                            if !(field in (:chord_uvec, :bl))
                                @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                            end
                        end

                        if positive_x_rotation
                            # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                            # If we're, applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
                            trans_phi = SteadyRotYTransformation(τ, 0.0, -twist)
                            chord_uvec_check = @SVector [0.0, 0.0, -1.0]
                        else
                            # If we're doing a negative-x rotation, we're applying the twist about the negative y axis.
                            # If we're applying the twist about the negative y axis, then we need to do a positive rotation about the y axis to undo it.
                            trans_phi = SteadyRotYTransformation(τ, 0.0, twist)
                            chord_uvec_check = @SVector [0.0, 0.0, 1.0]
                        end
                        se_no_twist = se |> trans_phi
                        @test se_no_twist.chord_uvec ≈ chord_uvec_check
                    end
                end
            end
        end
    end

end

end # module
