module CCBladeHelperTests

using AcousticAnalogies
using CCBlade
using DelimitedFiles
using KinematicCoordinateTransformations
using StaticArrays
using Test

include("gen_test_data/gen_ccblade_data/constants.jl")
using .CCBladeTestCaseConstants
ccbc = CCBladeTestCaseConstants

@testset "CCBlade private utils tests" begin

    Δr = (ccbc.Rtip - ccbc.Rhub)/10
    r = (ccbc.Rhub+0.5*Δr):Δr:(ccbc.Rtip-0.5*Δr)
    precone = 3*pi/180
    rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; precone=precone, turbine=false)

    dummy = similar(r)
    sections = Section.(r, dummy, dummy, dummy)

    @test all(AcousticAnalogies.get_ccblade_dradii(rotor, sections) .≈ Δr)
end

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
        Rhub = 0.2
        Rtip = 1.1
        num_blades = 3
        r = 2.0
        Δr = 0.1
        chord = 0.15
        twist = 15*pi/180
        af = nothing
        Vx = 35.0
        Vy = 100.0
        rho = 1.25
        Np = 2.0
        Tp = 3.0
        dummies = fill(0.0, 13)
        out = CCBlade.Outputs(Np, Tp, dummies...)
        section = CCBlade.Section(r, chord, twist, af)
        op = CCBlade.OperatingPoint(Vx, Vy, rho)
        area_per_chord2 = 0.1
        τ = 0.1
        rotor0precone = CCBlade.Rotor(Rhub, Rtip, num_blades; turbine=false, precone=0.0)
        for positive_x_rotation in [true, false]
            se_0theta0precone = CompactSourceElement(rotor0precone, section, op, out, 0.0, Δr, area_per_chord2, τ, positive_x_rotation)
            for precone in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                rotor = CCBlade.Rotor(Rhub, Rtip, num_blades; turbine=false, precone=precone)
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
                    se = TBLTESourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
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
        Rhub = 0.2
        Rtip = 1.1
        num_blades = 3
        r = 2.0
        Δr = 0.1
        chord = 0.15
        twist = 15*pi/180
        af = nothing
        Vx = 35.0
        Vy = 100.0
        rho = 1.25
        Np = 2.0
        Tp = 3.0
        dummies = fill(0.0, 13)
        out = CCBlade.Outputs(Np, Tp, dummies...)
        op = CCBlade.OperatingPoint(Vx, Vy, rho)
        τ = 0.1
        Δτ = 0.02
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        rotor0precone = CCBlade.Rotor(Rhub, Rtip, num_blades; turbine=false, precone=0.0)
        for positive_x_rotation in [true, false]
            for twist in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                section = CCBlade.Section(r, chord, twist, af)
                se_0theta0precone = TBLTESourceElement(rotor0precone, section, op, out, 0.0, Δr, τ, Δτ, bl, positive_x_rotation)
                for precone in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    rotor = CCBlade.Rotor(Rhub, Rtip, num_blades; turbine=false, precone=precone)
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

@testset "CCBlade CompactSourceElement test" begin

    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)
        # Read in the loading data.
        data = readdlm("gen_test_data/gen_ccblade_data/ccblade_omega11.csv", ',')
        Np = data[:, 1]
        Tp = data[:, 2]

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # Only care about getting the loading into the CCBlade Output structs.
        dummies = fill(0.0, 13)
        outs = Outputs.(Np, Tp, dummies...)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        aoc2 = fill(ccbc.area_over_chord_squared, length(sections))
        ses_helper = source_elements_ccblade(rotor, sections, ops, outs, aoc2, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way. First get the
        # transformation objects.
        rot_axis = @SVector [1.0, 0.0, 0.0]
        blade_axis = @SVector [0.0, 1.0, 0.0]
        y0_hub = @SVector [0.0, 0.0, 0.0]  # m
        v0_hub = ccbc.v.*rot_axis
        t0 = 0.0
        if positive_x_rotation
            rot_trans = SteadyRotXTransformation(t0, omega, 0.0)
        else
            rot_trans = SteadyRotXTransformation(t0, -omega, 0.0)
        end
        const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

        # Need the source times.
        dt = src_time_range/num_src_times
        src_times = t0 .+ (0:num_src_times-1).*dt

        # This is just an array of the angular offsets of each blade.
        θs = 2*pi/ccbc.num_blades.*(0:(ccbc.num_blades-1))

        # Radial spacing.
        dradii = get_dradii(ccbc.radii, ccbc.Rhub, ccbc.Rtip)

        # Reshape stuff for broadcasting.
        radii = reshape(ccbc.radii, 1, :, 1)
        dradii = reshape(dradii, 1, :, 1)
        cs_area = reshape(ccbc.area_over_chord_squared.*ccbc.chord.^2, 1, :, 1)
        Np = reshape(Np, 1, :, 1)
        Tp = reshape(Tp, 1, :, 1)
        src_times = reshape(src_times, :, 1, 1)  # This isn't really necessary.
        θs = reshape(θs, 1, 1, :)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Transform the source elements.
        if positive_x_rotation
            ses = CompactSourceElement.(ccbc.rho, ccbc.c0, radii, θs, dradii, cs_area, -Np, 0.0, Tp, src_times) .|> trans
        else
            ses = CompactSourceElement.(ccbc.rho, ccbc.c0, radii, θs, dradii, cs_area, -Np, 0.0, -Tp, src_times) .|> trans
        end

        for field in fieldnames(CompactSourceElement)
            @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
        end
    end
end

end
