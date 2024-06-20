module BroadbandSourceElementTests

using AcousticAnalogies
using AcousticAnalogies: calculate_bpm_test
using AcousticMetrics: AcousticMetrics
using CCBlade
using DelimitedFiles: DelimitedFiles
using FLOWMath: linear
using KinematicCoordinateTransformations
using StaticArrays
using JLD2
using Test
using LinearAlgebra: norm, dot, cross

include("gen_test_data/gen_ccblade_data/constants.jl")
using .CCBladeTestCaseConstants
ccbc = CCBladeTestCaseConstants

@testset "Twist and rotation tests" begin
    # So, the way this should work: first do the twist, then do the theta rotation.
    # The twist could be either about the positive y axis or negative y axis.
    # Then the theta rotation is always about the x axis.
    c0 = 1.1
    nu = 1.2
    r = 2.0
    Δr = 0.1
    chord = 1.3
    h = 1.4
    Psi = 1.5
    vn = 2.0
    vr = 3.0
    vc = 4.0
    τ = 0.1
    Δτ = 0.02
    bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
    blade_tip = 3.0 # should be a blade tip struct, but doesn't matter for these tests.
    for setype in [TBLTESourceElement, LBLVSSourceElement, TEBVSSourceElement, TipVortexSourceElement, CombinedNoTipBroadbandSourceElement, CombinedWithTipBroadbandSourceElement]
        for twist_about_positive_y in [true, false]
            if setype == CombinedWithTipBroadbandSourceElement
                se_0twist0theta = setype(c0, nu, r, 0.0, Δr, chord, 0.0, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
            elseif setype == TipVortexSourceElement
                se_0twist0theta = setype(c0, r, 0.0, Δr, chord, 0.0, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
            elseif setype in (TEBVSSourceElement, CombinedNoTipBroadbandSourceElement)
                se_0twist0theta = setype(c0, nu, r, 0.0, Δr, chord, 0.0, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
            else
                se_0twist0theta = setype(c0, nu, r, 0.0, Δr, chord, 0.0, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
            end
            for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

                for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    # The angle of attack depends on the twist and the fluid velocity
                    if twist_about_positive_y
                        alpha_check = ϕ - atan(-vn, -vc)
                    else
                        alpha_check = ϕ - atan(-vn, vc)
                    end
                    if setype == CombinedWithTipBroadbandSourceElement
                        se = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) |> trans_theta
                    elseif setype == TipVortexSourceElement
                        se = setype(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) |> trans_theta
                    elseif setype in (TEBVSSourceElement, CombinedNoTipBroadbandSourceElement)
                        se = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                    else
                        se = setype(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                    end
                    # Adjust the angles of attack to always be between -pi and pi.
                    alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                    alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                    @test alpha ≈ alpha_check

                    for field in fieldnames(setype)
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

                    # Make sure we get the same thing if we specify the velocity via a velocity magnitude and angle of attack.
                    # But need to make sure we use vr == 0.
                    if setype == CombinedWithTipBroadbandSourceElement
                        se_no_vr = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, 0.0, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
                    elseif setype == TipVortexSourceElement
                        se_no_vr = setype(c0, r, θ, Δr, chord, ϕ, vn, 0.0, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
                    elseif setype in (TEBVSSourceElement, CombinedNoTipBroadbandSourceElement)
                        se_no_vr = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, 0.0, vc, τ, Δτ, bl, twist_about_positive_y)
                    else
                        se_no_vr = setype(c0, nu, r, θ, Δr, chord, ϕ, vn, 0.0, vc, τ, Δτ, bl, twist_about_positive_y)
                    end
                    # Removing vr, the radial velocity component, shouldn't change the angle of attack.
                    alpha_no_vr = rem2pi(AcousticAnalogies.angle_of_attack(se_no_vr)+pi, RoundNearest) - pi
                    @test alpha_no_vr ≈ alpha_check
                    # Now create a source element using the velocity magnitude and angle of attack, check that we get the same thing.
                    U = sqrt(vn^2 + vc^2)
                    if setype == CombinedWithTipBroadbandSourceElement
                        se_from_U_α = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, U, alpha_no_vr, τ, Δτ, bl, blade_tip, twist_about_positive_y)
                    elseif setype == TipVortexSourceElement
                        se_from_U_α = setype(c0, r, θ, Δr, chord, ϕ, U, alpha_no_vr, τ, Δτ, bl, blade_tip, twist_about_positive_y)
                    elseif setype in (TEBVSSourceElement, CombinedNoTipBroadbandSourceElement)
                        se_from_U_α = setype(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, U, alpha_no_vr, τ, Δτ, bl, twist_about_positive_y)
                    else
                        se_from_U_α = setype(c0, nu, r, θ, Δr, chord, ϕ, U, alpha_no_vr, τ, Δτ, bl, twist_about_positive_y)
                    end
                    for field in fieldnames(setype)
                        @test getproperty(se_from_U_α, field) ≈ getproperty(se_no_vr, field)
                    end

                end
            end
        end
    end
end

@testset "TBLTESourceElement twist and rotation tests, CCBlade" begin
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
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CCBlade TBLTESourceElement complete test" begin
    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # What actually matters in the output structs are just W and phi.
        phi = range(45.0, 10.0; length=length(sections)) .* (pi/180)
        W = range(10.0, 11.0; length=length(sections))
        outs = Outputs.(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phi, 0.0, W, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        bls = [AcousticAnalogies.TrippedN0012BoundaryLayer()]
        ses_helper = tblte_source_elements_ccblade(rotor, sections, ops, outs, bls, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way.
        # First get the transformation objects.
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

        # Need the kinematic viscosity.
        nus = getproperty.(ops, :mu) ./ getproperty.(ops, :rho)

        # Also need the velocity in each direction.
        if positive_x_rotation
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. -W*cos(phi)
        else
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. W*cos(phi)
        end

        # Reshape stuff for broadcasting.
        radii_rs = reshape(ccbc.radii, 1, :, 1)
        dradii_rs = reshape(dradii, 1, :, 1)
        phi_rs = reshape(phi, 1, :, 1)
        W_rs = reshape(W, 1, :, 1)
        src_times_rs = reshape(src_times, :, 1, 1)  # This isn't really necessary.
        θs_rs = reshape(θs, 1, 1, :)
        nus_rs = reshape(nus, 1, :, 1)
        twist_rs = reshape(getproperty.(sections, :theta), 1, :, 1)
        chord_rs = reshape(getproperty.(sections, :chord), 1, :, 1)
        vn_rs = reshape(vn, 1, :, 1)
        vr_rs = reshape(vr, 1, :, 1)
        vc_rs = reshape(vc, 1, :, 1)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Transform the source elements.
        ses = TBLTESourceElement.(ccbc.c0, nus_rs, radii_rs, θs_rs, dradii_rs, chord_rs, twist_rs, vn_rs, vr_rs, vc_rs, src_times_rs, dt, bls, positive_x_rotation) .|> trans

        # Now check that we got the same thing.
        for field in fieldnames(TBLTESourceElement)
            if !(field in (:bl,))
                @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
            end
        end
    end
end

@testset "LBLVSSourceElement twist and rotation tests" begin
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
        se_0twist0theta = LBLVSSourceElement(c0, nu, r, 0.0, Δr, chord, 0.0, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

            for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                # The angle of attack depends on the twist and the fluid velocity
                if twist_about_positive_y
                    alpha_check = ϕ - atan(-vn, -vc)
                else
                    alpha_check = ϕ - atan(-vn, vc)
                end
                se = LBLVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                # Adjust the angles of attack to always be between -pi and pi.
                alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                @test alpha ≈ alpha_check

                for field in fieldnames(LBLVSSourceElement)
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

@testset "LBLVSSourceElement twist and rotation tests, CCBlade" begin
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
            se_0theta0precone = LBLVSSourceElement(rotor0precone, section, op, out, 0.0, Δr, τ, Δτ, bl, positive_x_rotation)
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
                    se = LBLVSSourceElement(rotor, section, op, out, θ, Δr, τ, Δτ, bl, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(LBLVSSourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if !(field in (:chord_uvec, :bl))
                            @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                        end
                    end

                    if positive_x_rotation
                        # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CCBlade LBLVSSourceElement complete test" begin
    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # What actually matters in the output structs are just W and phi.
        phi = range(45.0, 10.0; length=length(sections)) .* (pi/180)
        W = range(10.0, 11.0; length=length(sections))
        outs = Outputs.(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phi, 0.0, W, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        bls = [AcousticAnalogies.TrippedN0012BoundaryLayer()]
        ses_helper = lblvs_source_elements_ccblade(rotor, sections, ops, outs, bls, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way.
        # First get the transformation objects.
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

        # Need the kinematic viscosity.
        nus = getproperty.(ops, :mu) ./ getproperty.(ops, :rho)

        # Also need the velocity in each direction.
        if positive_x_rotation
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. -W*cos(phi)
        else
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. W*cos(phi)
        end

        # Reshape stuff for broadcasting.
        radii_rs = reshape(ccbc.radii, 1, :, 1)
        dradii_rs = reshape(dradii, 1, :, 1)
        phi_rs = reshape(phi, 1, :, 1)
        W_rs = reshape(W, 1, :, 1)
        src_times_rs = reshape(src_times, :, 1, 1)  # This isn't really necessary.
        θs_rs = reshape(θs, 1, 1, :)
        nus_rs = reshape(nus, 1, :, 1)
        twist_rs = reshape(getproperty.(sections, :theta), 1, :, 1)
        chord_rs = reshape(getproperty.(sections, :chord), 1, :, 1)
        vn_rs = reshape(vn, 1, :, 1)
        vr_rs = reshape(vr, 1, :, 1)
        vc_rs = reshape(vc, 1, :, 1)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Transform the source elements.
        ses = LBLVSSourceElement.(ccbc.c0, nus_rs, radii_rs, θs_rs, dradii_rs, chord_rs, twist_rs, vn_rs, vr_rs, vc_rs, src_times_rs, dt, bls, positive_x_rotation) .|> trans

        # Now check that we got the same thing.
        for field in fieldnames(LBLVSSourceElement)
            if !(field in (:bl,))
                @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
            end
        end
    end
end

@testset "TEBVSSourceElement twist and rotation tests" begin
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
    h = 0.1
    Psi = 0.2
    bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
    for twist_about_positive_y in [true, false]
        se_0twist0theta = TEBVSSourceElement(c0, nu, r, 0.0, Δr, chord, 0.0, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

            for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                # The angle of attack depends on the twist and the fluid velocity
                if twist_about_positive_y
                    alpha_check = ϕ - atan(-vn, -vc)
                else
                    alpha_check = ϕ - atan(-vn, vc)
                end
                se = TEBVSSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                # Adjust the angles of attack to always be between -pi and pi.
                alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                @test alpha ≈ alpha_check

                for field in fieldnames(TEBVSSourceElement)
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

@testset "TEBVSSourceElement twist and rotation tests, CCBlade" begin
    # Create the CCBlade objects.
    τ = 0.1
    Δτ = 0.02
    h = 0.1
    Psi = 0.2
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
            se_0theta0precone = TEBVSSourceElement(rotor0precone, section, op, out, 0.0, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation)
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
                    se = TEBVSSourceElement(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(TEBVSSourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if !(field in (:chord_uvec, :bl))
                            @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                        end
                    end

                    if positive_x_rotation
                        # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CCBlade TEBVSSourceElement complete test" begin
    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # What actually matters in the output structs are just W and phi.
        phi = range(45.0, 10.0; length=length(sections)) .* (pi/180)
        W = range(10.0, 11.0; length=length(sections))
        outs = Outputs.(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phi, 0.0, W, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        bls = [AcousticAnalogies.TrippedN0012BoundaryLayer()]
        hs = range(0.1, 0.2; length=length(sections))
        Psis = range(0.2, 0.3; length=length(sections))
        ses_helper = tebvs_source_elements_ccblade(rotor, sections, ops, outs, hs, Psis, bls, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way.
        # First get the transformation objects.
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

        # Need the kinematic viscosity.
        nus = getproperty.(ops, :mu) ./ getproperty.(ops, :rho)

        # Also need the velocity in each direction.
        if positive_x_rotation
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. -W*cos(phi)
        else
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. W*cos(phi)
        end

        # Reshape stuff for broadcasting.
        radii_rs = reshape(ccbc.radii, 1, :, 1)
        dradii_rs = reshape(dradii, 1, :, 1)
        phi_rs = reshape(phi, 1, :, 1)
        W_rs = reshape(W, 1, :, 1)
        src_times_rs = reshape(src_times, :, 1, 1)  # This isn't really necessary.
        θs_rs = reshape(θs, 1, 1, :)
        nus_rs = reshape(nus, 1, :, 1)
        twist_rs = reshape(getproperty.(sections, :theta), 1, :, 1)
        chord_rs = reshape(getproperty.(sections, :chord), 1, :, 1)
        vn_rs = reshape(vn, 1, :, 1)
        vr_rs = reshape(vr, 1, :, 1)
        vc_rs = reshape(vc, 1, :, 1)
        h_rs = reshape(hs, 1, :, 1)
        Psi_rs = reshape(Psis, 1, :, 1)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Transform the source elements.
        ses = TEBVSSourceElement.(ccbc.c0, nus_rs, radii_rs, θs_rs, dradii_rs, chord_rs, twist_rs, h_rs, Psi_rs, vn_rs, vr_rs, vc_rs, src_times_rs, dt, bls, positive_x_rotation) .|> trans

        # Now check that we got the same thing.
        for field in fieldnames(TEBVSSourceElement)
            if !(field in (:bl,))
                @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
            end
        end
    end
end

@testset "TipVortexSourceElement twist and rotation tests" begin
    # So, the way this should work: first do the twist, then do the theta rotation.
    # The twist could be either about the positive y axis or negative y axis.
    # Then the theta rotation is always about the x axis.
    c0 = 1.1
    # nu = 1.2
    r = 2.0
    Δr = 0.1
    chord = 1.3
    vn = 2.0
    vr = 3.0
    vc = 4.0
    τ = 0.1
    Δτ = 0.02
    bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
    blade_tip = 3.0 # should be a blade tip struct, but doesn't matter for these tests.
    for twist_about_positive_y in [true, false]
        se_0twist0theta = TipVortexSourceElement(c0, r, 0.0, Δr, chord, 0.0, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

            for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                # The angle of attack depends on the twist and the fluid velocity
                if twist_about_positive_y
                    alpha_check = ϕ - atan(-vn, -vc)
                else
                    alpha_check = ϕ - atan(-vn, vc)
                end
                se = TipVortexSourceElement(c0, r, θ, Δr, chord, ϕ, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) |> trans_theta
                # Adjust the angles of attack to always be between -pi and pi.
                alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                @test alpha ≈ alpha_check

                for field in fieldnames(TipVortexSourceElement)
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

@testset "TipVortexSourceElement twist and rotation tests, CCBlade" begin
    # Create the CCBlade objects.
    τ = 0.1
    Δτ = 0.02
    bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
    blade_tip = AcousticAnalogies.RoundedTip()
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
            se_0theta0precone = TipVortexSourceElement(rotor0precone, section, op, out, 0.0, Δr, τ, Δτ, bl, blade_tip, positive_x_rotation)
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
                    se = TipVortexSourceElement(rotor, section, op, out, θ, Δr, τ, Δτ, bl, blade_tip, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(TipVortexSourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if !(field in (:chord_uvec, :bl, :blade_tip))
                            @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                        end
                    end

                    if positive_x_rotation
                        # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CCBlade TipVortexSourceElement complete test" begin
    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # What actually matters in the output structs are just W and phi.
        phi = range(45.0, 10.0; length=length(sections)) .* (pi/180)
        W = range(10.0, 11.0; length=length(sections))
        outs = Outputs.(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phi, 0.0, W, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Radial spacing.
        dradii = get_dradii(ccbc.radii, ccbc.Rhub, ccbc.Rtip)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
        blade_tip = AcousticAnalogies.RoundedTip()
        ses_helper = tip_vortex_source_elements_ccblade(rotor, sections[end], ops[end], outs[end], dradii[end], bl, blade_tip, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way.
        # First get the transformation objects.
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

        # Need the kinematic viscosity.
        # nus = getproperty.(ops, :mu) ./ getproperty.(ops, :rho)

        # Also need the velocity in each direction.
        if positive_x_rotation
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. -W*cos(phi)
        else
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. W*cos(phi)
        end

        # Reshape stuff for broadcasting.
        θs_rs = reshape(θs, 1, 1, :)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Transform the source elements.
        ses = TipVortexSourceElement.(ccbc.c0, ccbc.radii[end], θs_rs, dradii[end], sections[end].chord, sections[end].theta, vn[end], vr[end], vc[end], src_times, dt, Ref(bl), Ref(blade_tip), positive_x_rotation) .|> trans

        # Now check that we got the same thing.
        for field in fieldnames(TipVortexSourceElement)
            if !(field in (:bl, :blade_tip))
                @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
            end
        end
    end
end

@testset "CombinedNoTipBroadbandSourceElement twist and rotation tests" begin
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
    h = 0.1
    Psi = 0.2
    bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
    for twist_about_positive_y in [true, false]
        se_0twist0theta = CombinedNoTipBroadbandSourceElement(c0, nu, r, 0.0, Δr, chord, 0.0, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

            for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                # The angle of attack depends on the twist and the fluid velocity
                if twist_about_positive_y
                    alpha_check = ϕ - atan(-vn, -vc)
                else
                    alpha_check = ϕ - atan(-vn, vc)
                end
                se = CombinedNoTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, twist_about_positive_y) |> trans_theta
                # Adjust the angles of attack to always be between -pi and pi.
                alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                @test alpha ≈ alpha_check

                for field in fieldnames(CombinedNoTipBroadbandSourceElement)
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

@testset "CombinedNoTipBroadbandSourceElement twist and rotation tests, CCBlade" begin
    # Create the CCBlade objects.
    τ = 0.1
    Δτ = 0.02
    h = 0.1
    Psi = 0.2
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
            se_0theta0precone = CombinedNoTipBroadbandSourceElement(rotor0precone, section, op, out, 0.0, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation)
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
                    se = CombinedNoTipBroadbandSourceElement(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(CombinedNoTipBroadbandSourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if !(field in (:chord_uvec, :bl))
                            @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                        end
                    end

                    if positive_x_rotation
                        # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CombinedWithTipBroadbandSourceElement twist and rotation tests" begin
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
    h = 0.1
    Psi = 0.2
    bl = 2.0 # should be a boundary layer struct, but doesn't matter for these tests.
    blade_tip = 3.0 # should be a blade tip struct, but doesn't matter for these tests.
    for twist_about_positive_y in [true, false]
        se_0twist0theta = CombinedWithTipBroadbandSourceElement(c0, nu, r, 0.0, Δr, chord, 0.0, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y)
        for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
            trans_theta = SteadyRotXTransformation(τ, 0.0, -θ)

            for ϕ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                # The angle of attack depends on the twist and the fluid velocity
                if twist_about_positive_y
                    alpha_check = ϕ - atan(-vn, -vc)
                else
                    alpha_check = ϕ - atan(-vn, vc)
                end
                se = CombinedWithTipBroadbandSourceElement(c0, nu, r, θ, Δr, chord, ϕ, h, Psi, vn, vr, vc, τ, Δτ, bl, blade_tip, twist_about_positive_y) |> trans_theta
                # Adjust the angles of attack to always be between -pi and pi.
                alpha_check = rem2pi(alpha_check+pi, RoundNearest) - pi
                alpha = rem2pi(AcousticAnalogies.angle_of_attack(se)+pi, RoundNearest) - pi
                @test alpha ≈ alpha_check

                for field in fieldnames(CombinedWithTipBroadbandSourceElement)
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

@testset "CombinedWithTipBroadbandSourceElement twist and rotation tests, CCBlade" begin
    # Create the CCBlade objects.
    τ = 0.1
    Δτ = 0.02
    h = 0.1
    Psi = 0.2
    bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
    blade_tip = AcousticAnalogies.RoundedTip()
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
            se_0theta0precone = CombinedWithTipBroadbandSourceElement(rotor0precone, section, op, out, 0.0, Δr, h, Psi, τ, Δτ, bl, blade_tip, positive_x_rotation)
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
                    se = CombinedWithTipBroadbandSourceElement(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, blade_tip, positive_x_rotation) |> trans
                    # Check that we got the same thing:
                    for field in fieldnames(CombinedWithTipBroadbandSourceElement)
                        # The twist changes the unit vector in the chord direction, but nothing else, so ignore that for now.
                        if !(field in (:chord_uvec, :bl, :blade_tip))
                            @test getproperty(se, field) ≈ getproperty(se_0theta0precone, field)
                        end
                    end

                    if positive_x_rotation
                        # If we're doing a positive-x rotation, we're applying the twist about the positive y axis.
                        # If we're applying the twist about the positive y axis, then we need to do a negative rotation about the y axis to undo it.
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

@testset "CCBlade combined broadband source elements complete test" begin
    for positive_x_rotation in [true, false]
        omega = 2200*(2*pi/60)

        # Create the CCBlade objects.
        rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
        sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
        ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
        # What actually matters in the output structs are just W and phi.
        phi = range(45.0, 10.0; length=length(sections)) .* (pi/180)
        W = range(10.0, 11.0; length=length(sections))
        outs = Outputs.(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, phi, 0.0, W, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Set the source time stuff.
        num_blade_passes = 3
        steps_per_blade_pass = 8
        num_src_times = num_blade_passes*steps_per_blade_pass
        bpp = 2*pi/omega/ccbc.num_blades
        src_time_range = num_blade_passes*bpp

        # Finally get all the source elements.
        # bls = [AcousticAnalogies.TrippedN0012BoundaryLayer()]
        # bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
        # bls = Fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), num_radial)
        bls = fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), length(sections))
        hs = range(0.1, 0.2; length=length(sections))
        Psis = range(0.2, 0.3; length=length(sections))
        blade_tip = AcousticAnalogies.RoundedTip()
        ses_no_tip_helper, ses_with_tip_helper = combined_broadband_source_elements_ccblade(rotor, sections, ops, outs, hs, Psis, bls, blade_tip, src_time_range, num_src_times, positive_x_rotation)

        # Now need to get the source elements the "normal" way.
        # First get the transformation objects.
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

        # Need the kinematic viscosity.
        nus = getproperty.(ops, :mu) ./ getproperty.(ops, :rho)

        # Also need the velocity in each direction.
        if positive_x_rotation
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. -W*cos(phi)
        else
            vn = @. -W*sin(phi)
            vr = zeros(eltype(vn), length(vn))
            vc = @. W*cos(phi)
        end

        # Reshape stuff for broadcasting.
        radii_rs = reshape(ccbc.radii, 1, :, 1)
        dradii_rs = reshape(dradii, 1, :, 1)
        phi_rs = reshape(phi, 1, :, 1)
        W_rs = reshape(W, 1, :, 1)
        # src_times_rs = reshape(src_times, :, 1, 1)  # This isn't really necessary.
        θs_rs = reshape(θs, 1, 1, :)
        nus_rs = reshape(nus, 1, :, 1)
        twist_rs = reshape(getproperty.(sections, :theta), 1, :, 1)
        chord_rs = reshape(getproperty.(sections, :chord), 1, :, 1)
        hs_rs = reshape(hs, 1, :, 1)
        Psis_rs = reshape(Psis, 1, :, 1)
        vn_rs = reshape(vn, 1, :, 1)
        vr_rs = reshape(vr, 1, :, 1)
        vc_rs = reshape(vc, 1, :, 1)
        bls_rs = reshape(bls, 1, :, 1)

        # Get all the transformations.
        trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

        # Now need to split things into the with tip and no tip stuff.
        radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
        dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
        phi_rs_no_tip = @view phi_rs[:, begin:end-1, :]
        W_rs_no_tip = @view W_rs[:, begin:end-1, :]
        nus_rs_no_tip = @view nus_rs[:, begin:end-1, :]
        twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
        chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
        hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
        Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
        vn_rs_no_tip = @view vn_rs[:, begin:end-1, :]
        vr_rs_no_tip = @view vr_rs[:, begin:end-1, :]
        vc_rs_no_tip = @view vc_rs[:, begin:end-1, :]
        bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

        radii_rs_with_tip = @view radii_rs[:, end:end, :]
        dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
        phi_rs_with_tip = @view phi_rs[:, end:end, :]
        W_rs_with_tip = @view W_rs[:, end:end, :]
        nus_rs_with_tip = @view nus_rs[:, end:end, :]
        twist_rs_with_tip = @view twist_rs[:, end:end, :]
        chord_rs_with_tip = @view chord_rs[:, end:end, :]
        hs_rs_with_tip = @view hs_rs[:, end:end, :]
        Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
        vn_rs_with_tip = @view vn_rs[:, end:end, :]
        vr_rs_with_tip = @view vr_rs[:, end:end, :]
        vc_rs_with_tip = @view vc_rs[:, end:end, :]
        bls_rs_with_tip = @view bls_rs[:, end:end, :]

        # Transform the source elements.
        ses_no_tip = CombinedNoTipBroadbandSourceElement.(ccbc.c0, nus_rs_no_tip, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord_rs_no_tip, twist_rs_no_tip, hs_rs_no_tip, Psis_rs_no_tip, vn_rs_no_tip, vr_rs_no_tip, vc_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans

        ses_with_tip = CombinedWithTipBroadbandSourceElement.(ccbc.c0, nus_rs_with_tip, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord_rs_with_tip, twist_rs_with_tip, hs_rs_with_tip, Psis_rs_with_tip, vn_rs_with_tip, vr_rs_with_tip, vc_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

        # Now check that we got the same thing.
        for field in fieldnames(CombinedNoTipBroadbandSourceElement)
            if !(field in (:bl,))
                @test all(getproperty.(ses_no_tip_helper, field) .≈ getproperty.(ses_no_tip, field))
            end
        end
        for field in fieldnames(CombinedWithTipBroadbandSourceElement)
            if !(field in (:bl, :blade_tip))
                @test all(getproperty.(ses_with_tip_helper, field) .≈ getproperty.(ses_with_tip, field))
            end
        end
    end
end

@testset "directivity function tests" begin
    # None of this stuff matters for the directivity functions.
    c0 = 2.0
    nu = 3.0
    dr = 0.1
    chord = 1.1
    τ = 0.2
    dτ = 0.01
    bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
    chord_cross_span_to_get_top_uvec = true

    # This stuff actually matters.
    # Need the fluid velocity to be zero so we can ignore the denominator.
    y1dot = @SVector [0.0, 0.0, 0.0]
    y1dot_fluid = @SVector [0.0, 0.0, 0.0]
    y0dot = @SVector [0.0, 0.0, 0.0]
    chord_uvec = [1.0, 0.0, 0.0]
    span_uvec = [0.0, 1.0, 0.0]

    # Create a source element.
    se = AcousticAnalogies.TBLTESourceElement(c0, nu, dr, chord, y0dot, y1dot, y1dot_fluid, τ, dτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)

    # Now, create an observer at different places, and check that we get the correct directivity function.
    for x_er in [-5.0, -3.0, 1.5, 4.0]
        for y_er in [-5.0, -3.0, 1.5, 4.0]
            for z_er in [-5.0, -3.0, 1.5, 4.0]
                x = @SVector [x_er, y_er, z_er]
                obs = AcousticAnalogies.StationaryAcousticObserver(x)
                r_er_check = sqrt(x_er^2 + y_er^2 + z_er^2)
                Θ_er = acos(x_er/r_er_check)
                Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2))

                # Observer time doesn't matter since the observer is stationary.
                t_obs = 7.0
                x_obs = obs(t_obs)
                Dl_check = (sin(Θ_er)^2) * (sin(Φ_er)^2)
                Dh_check = 2*(sin(0.5*Θ_er)^2) * (sin(Φ_er)^2)
                top_is_suction = true
                r_er, Dl, Dh = AcousticAnalogies.directivity(se, obs(t_obs), top_is_suction)
                @test r_er ≈ r_er_check
                @test Dl ≈ Dl_check
                @test Dh ≈ Dh_check

                # Now, rotate and translate both the source and the observer.
                # The directivity functions should be the same.
                # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
                trans1 = SteadyRotXTransformation(t_obs, 0.0, 3.0*pi/180)
                trans2 = SteadyRotYTransformation(t_obs, 0.0, 4.0*pi/180)
                trans3 = SteadyRotZTransformation(t_obs, 0.0, 5.0*pi/180)
                # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
                x_trans = @SVector [2.0, 3.0, 4.0]
                v_trans = @SVector [0.0, 0.0, 0.0]
                trans4 = ConstantVelocityTransformation(t_obs, x_trans, v_trans)

                # Transform the source and observer.
                trans = compose(t_obs, trans4, compose(t_obs, trans3, compose(t_obs, trans2, trans1)))
                se_trans = trans(se)
                obs_trans = AcousticAnalogies.StationaryAcousticObserver(trans(t_obs, obs(t_obs)))

                # Check that the directivity functions didn't change.
                r_er, Dl, Dh = AcousticAnalogies.directivity(se_trans, obs_trans(t_obs), top_is_suction)
                @test r_er ≈ r_er_check
                @test Dl ≈ Dl_check
                @test Dh ≈ Dh_check
            end
        end
    end

end

@testset "angle of attack test" begin

    @testset "TBLTESourceElement" begin
        for twist_about_positive_y in [true, false]
            # None of this stuff matters for the angle of attack.
            c0 = 2.0
            nu = 3.0
            dr = 0.1
            chord = 1.1
            dτ = 0.01
            r = 0.5
            omega = 101.0
            bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

            # We'll create a random transformation to check that rotating and displacing the source doesn't change the angle of attack.
            # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
            τ = 0.2
            trans1 = SteadyRotXTransformation(τ, 0.0, 3.0*pi/180)
            trans2 = SteadyRotYTransformation(τ, 0.0, 4.0*pi/180)
            trans3 = SteadyRotZTransformation(τ, 0.0, 5.0*pi/180)
            # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
            x_trans = @SVector [2.0, 3.0, 4.0]
            v_trans = @SVector [0.0, 0.0, 0.0]
            trans4 = ConstantVelocityTransformation(τ, x_trans, v_trans)
            # Combine all the transformations into one.
            trans = compose(τ, trans4, compose(τ, trans3, compose(τ, trans2, trans1)))

            # This stuff does matter for angle of attack.
            Vx = 5.5
            u = 0.1*Vx
            Vy = omega*r
            v = 0.05*Vy
            # So, let's say I'm in the usual frame of reference: moving axially in the positive x direction, rotating about the positive x axis if `twist_about_positive_y` is `true`, negative x axis if `twist_about_positive_y` is `false`, initially aligned with the y axis.
            # Then, from the perspective of the blade element, the fluid velocity in the axial direction (`x`) is `-(Vx + u)`, and the velocity in the tangential direction is `(-Vy + v)`.
            #
            # vr shouldn't matter at all, so check that.
            for use_induction in [true, false]
                for vr in [0.0, 2.1, -2.1]
                    for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                        for twist in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                            for vn_sign in [1, -1]
                                for vc_sign in [1, -1]
                                    vn = -(Vx + u)*vn_sign
                                    vc = (-Vy + v)*vc_sign
                                    se = AcousticAnalogies.TBLTESourceElement{AcousticAnalogies.BrooksBurleyDirectivity,use_induction,AcousticAnalogies.NoMachCorrection,true}(c0, nu, r, θ, dr, chord, twist, vn, vr, vc, τ, dτ, bl, twist_about_positive_y)
                                    # And then the angle of attack will be `twist - atan(-vn, -vc)`, where `twist` is the twist of the blade, `vn` is the velocity in the axial direction, and `vc` is the velocity in the circumferential/tangential direction.
                                    # But we need to switch the direction of the velocity vector, since I'm thinking of them in opposite directions (eg the angle of attack is zero when the velocity and chordwise vector from trailing edge to leading edge are exactly opposite).
                                    # The rem2pi will give us back an equivalent angle in the interval [-π, π].
                                    if twist_about_positive_y
                                        # If the twist is about the positive y axis, then the assumption is that the twist and velocity angles are zero when they are aligned with positive z axis.
                                        # In the "usual" operation, the axial component of the blade-fixed frame velocity will be in the negative x direction, and the circumferential component of the blade-fixed frame velocity will be in the negative z direction.
                                        # So we need to switch both of those signs to get the correct angle.
                                        angle_of_attack_check = rem2pi(twist - atan(-vn, -vc), RoundNearest)
                                    else
                                        # If the twist is about the negative y axis, then the assumption is the twist and velocity angles are zero when they are aligned with the negative z axis.
                                        # In then "usual" operation, the axial component of the blade-fixed frame velocity will be in the negative x direction, and the circumferential component of the blade fixed frame velocity will be in the positive z direction.
                                        # So we need to switch the sign of the axial component, but not the circumferential.
                                        angle_of_attack_check = rem2pi(twist - atan(-vn, vc), RoundNearest)
                                    end

                                    @test AcousticAnalogies.angle_of_attack(se) ≈ angle_of_attack_check

                                    if use_induction
                                        # Flow speed normal to span, including induction:
                                        U_check = sqrt(vn^2 + vc^2)
                                    else
                                        # Flow speed normal to span, not including induction:
                                        U_check = sqrt(Vx^2 + Vy^2)
                                    end

                                    if use_induction
                                        # If we're including induction in the flow speed normal to span calculation, we can check it now.
                                        @test AcousticAnalogies.speed_normal_to_span(se) ≈ U_check
                                    end

                                    # Apply the random transformation we created.
                                    se_trans = trans(se)

                                    # Make sure we get the same angle of attack as before.
                                    @test AcousticAnalogies.angle_of_attack(se_trans) ≈ angle_of_attack_check

                                    # Now, instead of doing a transformation that just displaces the source element, let's do one that changes the velocity.
                                    # So, how are we going to transform the source element into the fluid frame?
                                    # Well, first we say we're rotating around the x axix at a rate ω or -ω, depending on the value of `twist_about_positive_y`
                                    # Oh, but wait.
                                    # We're switching the sign.
                                    # Hmm... what to do about that?
                                    # Well, the definition of the fluid frame is one that has the "freestream velocity" as zero.
                                    # So that, I think, means we need to remove the axial and circumferential (rotational) velocity.
                                    # So remove Vx and Vy.
                                    # I should be able to do that.
                                    # So, first, let's think about getting rid of Vy.
                                    # If we rotate about the positive x axis, then that will increase the velocity in the z direction (since the blade is initially aligned with the y axis).
                                    # If we rotate about the negative x axis, then that will decrease the velocity in the z direction (again, since the blade is initially aligned with the y axis).
                                    # OK, then.
                                    trans_rot = KinematicCoordinateTransformations.SteadyRotXTransformation(τ, omega*vc_sign, 0.0)

                                    # Now, for the x velocity, we just want to remove the Vx.
                                    x0 = @SVector [0.0, 0.0, 0.0]
                                    v0 = @SVector [Vx*vn_sign, 0.0, 0.0]
                                    trans_freestream = KinematicCoordinateTransformations.ConstantVelocityTransformation(τ, x0, v0)

                                    # Now compose the two transformations, and apply them to the source element.
                                    trans_global = compose(τ, trans_freestream, trans_rot)
                                    se_global = trans_global(se)

                                    # The angle of attack should be the same.
                                    @test AcousticAnalogies.angle_of_attack(se_global) ≈ angle_of_attack_check

                                    # Also, we should know what the source element and fluid velocities are, right?
                                    y1dot_check = @SVector [Vx*vn_sign, Vy*vc_sign*(-sin(θ)), Vy*vc_sign*(cos(θ))]
                                    y1dot_fluid_check = @SVector [-u*vn_sign, vr*cos(θ) + v*vc_sign*(-sin(θ)), vr*sin(θ) + v*vc_sign*(cos(θ))]

                                    @test se_global.y1dot ≈ y1dot_check
                                    @test se_global.y1dot_fluid ≈ y1dot_fluid_check

                                    # The flow speed normal to span, including induction or not, shouldn't have changed, so check that.
                                    @test AcousticAnalogies.speed_normal_to_span(se_global) ≈ U_check
                                end
                            end
                        end
                    end
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
        out, section, Δr, op, rotor = nothing, nothing, nothing, nothing, nothing
        jldopen(ccblade_fname, "r") do f
            out = f["outs"][1]
            section = f["sections"][1]
            Δr = f["sections"][2].r - f["sections"][1].r
            op = f["ops"][1]
            rotor = f["rotor"]
            @test rotor.precone ≈ 0.0
            @test op.pitch ≈ 0.0
        end
        for positive_x_rotation in [true, false]
            for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                se = AcousticAnalogies.TBLTESourceElement(rotor, section, op, out, θ, Δr, τ, Δτ, bl, positive_x_rotation)
                # The `chord_uvec` vector points from leading edge to trailing edge.
                # So we should be able to use that to figure out the angle it makes with the tangential/rotation direction.
                # The tangential/rotation direction can be found by crossing the the rotation axis with the position vector.
                # Then I can dot the chord with that direction, and the forward velocity axis, then use the arctan function to get the angle.
                if positive_x_rotation
                    rot_axis = @SVector [1, 0, 0]
                else
                    rot_axis = @SVector [-1, 0, 0]
                end
                tan_axis_tmp = cross(rot_axis, se.y0dot)
                tan_axis = tan_axis_tmp / norm(tan_axis_tmp)
                forward_axis = @SVector [1, 0, 0]
                # I'm visualizing the chord vector as going from trailing edge to leading edge, but it's leading edge to trailing edge in the TBLTESourceElement struct, so switch that.
                te_to_le = -se.chord_uvec
                twist_check = atan(dot(te_to_le, forward_axis), dot(te_to_le, tan_axis))
                @test twist_check ≈ section.theta

                # The angle of attack that AcousticAnalogies.jl calculates should match what CCBlade.jl has.
                @test AcousticAnalogies.angle_of_attack(se) ≈ out.alpha

                # Now, rotate and translate the source, which shouldn't change the twist or angle of attack, as long as we don't do anything that would change the velocity.
                # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
                trans1 = SteadyRotXTransformation(τ, 0.0, 3.0*pi/180)
                trans2 = SteadyRotYTransformation(τ, 0.0, 4.0*pi/180)
                trans3 = SteadyRotZTransformation(τ, 0.0, 5.0*pi/180)
                # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
                x_trans = @SVector [2.0, 3.0, 4.0]
                v_trans = @SVector [0.0, 0.0, 0.0]
                trans4 = ConstantVelocityTransformation(τ, x_trans, v_trans)

                # Transform the source.
                trans = compose(τ, trans4, compose(τ, trans3, compose(τ, trans2, trans1)))
                se_trans = trans(se)

                # Angle of attack should still be the same.
                @test AcousticAnalogies.angle_of_attack(se_trans) ≈ out.alpha

                # It'd be nice to check the twist too, but if that was wrong, the angle of attack would be wrong too.
                # And there are other tests already for `chord_uvec`.

                # Could I put the source element in the "fluid/global" frame now?
                # I'd need to know the forward velocity and rotation rate.
                # Well, the forward velocity would be Vx.
                Vx = op.Vx
                # And the rotation rate would be Vy/r.
                omega = op.Vy / section.r
                # Now, need to undo the rotation, which depends on `positive_x_rotation`.
                if positive_x_rotation
                    # If we're rotating about the positive x axis, need to apply a rotation around the negative x axis to undo it.
                    trans_rot = KinematicCoordinateTransformations.SteadyRotXTransformation(τ, -omega, 0.0)
                else
                    # If we're rotating about the negative x axis, need to apply a rotation around the positive x axis to undo it.
                    trans_rot = KinematicCoordinateTransformations.SteadyRotXTransformation(τ, omega, 0.0)
                end
                # Now, I'm assuming that the freestream/axial velocity is in the -x direction, so to undo that, move it in the positive x direction.
                x0 = @SVector [0.0, 0.0, 0.0]
                v0 = @SVector [Vx, 0.0, 0.0]
                trans_freestream = KinematicCoordinateTransformations.ConstantVelocityTransformation(τ, x0, v0)

                # Now compose the two transformations, and apply them to the original source element.
                trans_global = compose(τ, trans_freestream, trans_rot)
                se_global = trans_global(se)

                # The angle of attack should be the same.
                @test AcousticAnalogies.angle_of_attack(se_global) ≈ out.alpha
            end
        end
    end
end

@testset "BPM Report tests" begin

    @testset "BPM Figure 11a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 30.48e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        M_c = 0.8*M
        alphastar = 0.0
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure11-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1] # This is in kHz.
        SPL_s = bpm[:, 2]

        # At zero angle of attack the pressure and suction side predictions are the same.
        f_p = f_s
        SPL_p = SPL_s

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.029

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.029

                # These should all be very negative, since alphastar is zero:
                @test all(SPL_alpha_jl .< -100)
            end
        end
    end

    @testset "BPM Figure 11d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 30.48e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure11-d-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1] # This is in kHz.
        SPL_s = bpm[:, 2]

        # At zero angle of attack the pressure and suction side predictions are the same.
        f_p = f_s
        SPL_p = SPL_s

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.015

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.015

                # These should all be very negative, since alphastar is zero:
                @test all(SPL_alpha_jl .< -100)
            end
        end
    end

    @testset "BPM Figure 12a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 30.48e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 1.5*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure12-U71.3-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure12-U71.3-TBL-TE-pressure.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure12-U71.3-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.022

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.017

                SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                vmin, vmax = extrema(SPL_alpha)
                err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                @test maximum(err) < 0.037
            end
        end
    end

    @testset "BPM Figure 26a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure26-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # At zero angle of attack the pressure and suction side predictions are the same.
        f_p = f_s
        SPL_p = SPL_s

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.015

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.015

                # These should all be very negative, since alphastar is zero:
                @test all(SPL_alpha_jl .< -100)
            end
        end
    end

    @testset "BPM Figure 26d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure26-d-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # At zero angle of attack the pressure and suction side predictions are the same.
        f_p = f_s
        SPL_p = SPL_s

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.032

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.032

                # These should all be very negative, since alphastar is zero:
                @test all(SPL_alpha_jl .< -100)
            end
        end
    end

    @testset "BPM Figure 28a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 6.7*pi/180
        # Using the tripped boundary layer in this case.
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-a-TBL-TE-pressure.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-a-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.036

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.075

                SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                vmin, vmax = extrema(SPL_alpha)
                err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                @test maximum(err) < 0.039
            end
        end
    end

    @testset "BPM Figure 28d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 6.7*pi/180
        # Using the tripped boundary layer in this case.
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-d-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-d-TBL-TE-pressure.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure28-d-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.021

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.042

                SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                vmin, vmax = extrema(SPL_alpha)
                err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                @test maximum(err) < 0.040
            end
        end
    end

    @testset "BPM Figure 38d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 2.54e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        # Using the tripped boundary layer in this case.
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure38-d-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # At zero angle of attack the pressure and suction side predictions are the same.
        f_p = f_s
        SPL_p = SPL_s

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.026

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.026

                # These should all be very negative, since alphastar is zero:
                @test all(SPL_alpha_jl .< -100)
            end
        end
    end


    @testset "BPM Figure 39d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 2.54e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 4.8*pi/180
        # Using the tripped boundary layer in this case.
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure39-d-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure39-d-TBL-TE-pressure.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure39-d-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.036

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.043

                SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                vmin, vmax = extrema(SPL_alpha)
                err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                @test maximum(err) < 0.039
            end
        end
    end

    @testset "BPM Figure 45a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 30.48e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 1.5*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure45-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure45-a-TBL-TE-pressure.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure45-a-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure45-a-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                # # Now compare...
                # SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                # vmin, vmax = extrema(SPL_s)
                # err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                # @test maximum(err) < 0.036

                # SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                # vmin, vmax = extrema(SPL_p)
                # err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                # @test maximum(err) < 0.043

                # SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                # vmin, vmax = extrema(SPL_alpha)
                # err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                # @test maximum(err) < 0.039

                # The agreement with these ones aren't so great.
                    # Might be better if I grabbed the listing in the BPM appendix?
                # Might be better if I grabbed the listing in the BPM appendix?
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.037

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.058

                SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
                vmin, vmax = extrema(SPL_alpha)
                err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
                @test maximum(err) < 0.091

                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.053
            end
        end
    end

    @testset "BPM Figure 48c" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 22.86e-2  # chord in meters
        U = 39.6  # freestream velocity in m/s
        M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure48-c-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.083
            end
        end
    end

    @testset "BPM Figure 54a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 15.24e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 2.7*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure54-a-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.026
            end
        end
    end

    @testset "BPM Figure 59c" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 39.6  # freestream velocity in m/s
        M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure59-c-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.11
            end
        end
    end

    @testset "BPM Figure 60c" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 39.6  # freestream velocity in m/s
        M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 3.3*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure60-c-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.12
            end
        end
    end

    @testset "BPM Figure 60d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 10.16e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 3.3*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure60-d-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.026
            end
        end
    end

    @testset "BPM Figure 65d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 5.08e-2  # chord in meters
        U = 31.7  # freestream velocity in m/s
        M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure65-d-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test maximum(err) < 0.021
            end
        end
    end

    @testset "BPM Figure 66b" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 5.08e-2  # chord in meters
        U = 39.6  # freestream velocity in m/s
        M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 4.2*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure66-b-LBL-VS.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_lbl_vs = bpm[:, 1]
        SPL_lbl_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_lblvs=true, use_Ualpha=use_Ualpha)

                SPL_lbl_vs_jl_interp = linear(freqs, SPL_lbl_vs_jl, f_lbl_vs.*1e3)
                vmin, vmax = extrema(SPL_lbl_vs)
                err = abs.(SPL_lbl_vs_jl_interp .- SPL_lbl_vs)./(vmax - vmin)
                @test length(err) == 3
                @test err[1] < 0.089
                @test err[2] < 0.373
                @test err[3] < 0.746
            end
        end
    end

    @testset "BPM Figure 69a" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 5.08e-2  # chord in meters
        U = 71.3  # freestream velocity in m/s
        M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 15.4*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure69-a-separation.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_alpha = bpm[:, 1]
        SPL_alpha = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl)

            # Now compare...
            @test all(SPL_s_jl .≈ -100)
            @test all(SPL_p_jl .≈ -100)

            SPL_alpha_jl_interp = linear(freqs, SPL_alpha_jl, f_alpha.*1e3)
            vmin, vmax = extrema(SPL_alpha)
            err = abs.(SPL_alpha_jl_interp .- SPL_alpha)./(vmax - vmin)
            @test maximum(err) < 0.033
        end
    end

    @testset "BPM Figure 91" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 30.48e-2  # span in meters
        chord = 15.24e-2  # chord in meters
        speedofsound = 340.46
        U = 71.3  # freestream velocity in m/s
        # M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
        M = U/speedofsound
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        # alphatip = 0.71*10.8*pi/180
        alphastar = 10.8*pi/180
        bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
        # blade_tip = AcousticAnalogies.RoundedTip{AcousticAnalogies.BPMTipAlphaCorrection}()
        blade_tip = AcousticAnalogies.RoundedTip(AcousticAnalogies.BPMTipAlphaCorrection(), 0.0)

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure91-tip.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_tip = bpm[:, 1]
        SPL_tip = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_tip_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tip_vortex=true, blade_tip=blade_tip, use_Ualpha=use_Ualpha)

                SPL_tip_jl_interp = linear(freqs, SPL_tip_jl, f_tip.*1e3)
                vmin, vmax = extrema(SPL_tip)
                err = abs.(SPL_tip_jl_interp .- SPL_tip)./(vmax - vmin)
                @test maximum(err) < 0.047
            end
        end
    end

    @testset "BPM Figure 98b" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 69.5  # freestream velocity in m/s
        M = U/340.46
        h = 1.1e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        # Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
        # And TBL-TE is the only significant source, other than bluntness.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-b-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-2]) < 0.052
                @test maximum(err[1:end-2]) < 0.052
                @test maximum(err[1:end-1]) < 0.060
                @test maximum(err) < 0.171
            end
        end
    end

    @testset "BPM Figure 98c" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 69.5  # freestream velocity in m/s
        M = U/340.46
        h = 1.9e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        M_c = 0.8*M
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        # Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
        # And TBL-TE is the only significant source, other than bluntness.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-c-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-2]) < 0.040
                @test maximum(err[1:end-1]) < 0.189
                @test err[end] < 0.111
            end
        end
    end

    @testset "BPM Figure 98d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 69.5  # freestream velocity in m/s
        M = U/340.46
        h = 2.5e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        # Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
        # And TBL-TE is the only significant source, other than bluntness.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-a-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure98-d-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.053

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-2]) < 0.044
                @test err[end-1] < 0.089
                @test err[end] < 0.089
            end
        end
    end

    @testset "BPM Figure 99b" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 38.6  # freestream velocity in m/s
        M = U/340.46
        h = 1.1e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-2]) < 0.091
                @test         err[  end-1]  < 0.251
                @test         err[  end  ]  < 0.400
            end
        end
    end

    @testset "BPM Figure 99c" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 38.6  # freestream velocity in m/s
        M = U/340.46
        h = 1.9e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        # Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
        # And TBL-TE is the only significant source, other than bluntness.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-c-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-2]) < 0.057
                @test         err[  end-1]  < 0.070
                @test         err[  end  ]  < 0.256
            end
        end
    end

    @testset "BPM Figure 99d" begin
        nu = 1.4529e-5  # kinematic viscosity, m^2/s
        L = 45.72e-2  # span in meters
        chord = 60.96e-2  # chord in meters
        U = 38.6  # freestream velocity in m/s
        M = U/340.46
        h = 2.5e-3  # trailing edge bluntness in meters
        Psi = 14*pi/180  # bluntness angle in radians
        r_e = 1.22 # radiation distance in meters
        θ_e = 90*pi/180 
        Φ_e = 90*pi/180
        alphastar = 0.0*pi/180
        bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

        # Now, need to get the data from the BPM report.
        # Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
        # And TBL-TE is the only significant source, other than bluntness.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_s = bpm[:, 1]
        SPL_s = bpm[:, 2]

        # Suction and pressure are the same for zero angle of attack.
        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-b-TBL-TE-suction.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_p = bpm[:, 1]
        SPL_p = bpm[:, 2]

        fname = joinpath(@__DIR__, "bpm_data", "19890016302-figure99-d-bluntness.csv")
        bpm = DelimitedFiles.readdlm(fname, ',')
        f_teb_vs = bpm[:, 1]
        SPL_teb_vs = bpm[:, 2]

        for angle_of_attack_sign in [1, -1]
            for use_Ualpha in [false, true]
                freqs, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, angle_of_attack_sign*alphastar, bl; do_tebvs=true, h=h, Psi=Psi, use_Ualpha=use_Ualpha)

                # Now compare...
                SPL_s_jl_interp = linear(freqs, SPL_s_jl, f_s.*1e3)
                vmin, vmax = extrema(SPL_s)
                err = abs.(SPL_s_jl_interp .- SPL_s)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_p_jl_interp = linear(freqs, SPL_p_jl, f_p.*1e3)
                vmin, vmax = extrema(SPL_p)
                err = abs.(SPL_p_jl_interp .- SPL_p)./(vmax - vmin)
                @test maximum(err) < 0.077

                SPL_teb_vs_jl_interp = linear(freqs, SPL_teb_vs_jl, f_teb_vs.*1e3)
                vmin, vmax = extrema(SPL_teb_vs)
                err = abs.(SPL_teb_vs_jl_interp .- SPL_teb_vs)./(vmax - vmin)
                # Last two points are off.
                # Not sure why.
                @test maximum(err[1:end-3]) < 0.047
                @test err[end-2] < 0.068
                @test err[end-1] < 0.213
                @test err[end] < 0.225
            end
        end
    end
end

end # module
