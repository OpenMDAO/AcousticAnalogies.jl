module TBLTESourceElementTests

using AcousticAnalogies
using CCBlade
# using DelimitedFiles
using KinematicCoordinateTransformations
using StaticArrays
using JLD2
using Test
using LinearAlgebra: norm, dot, cross

@testset "TBLTESourceElement twist and rotation tests" begin
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

@testset "directivity function tests" begin
    # None of this stuff matters for the directivity functions.
    c0 = 2.0
    nu = 3.0
    dr = 0.1
    chord = 1.1
    τ = 0.2
    dτ = 0.01
    bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
    chord_cross_span_to_get_suction_uvec = true

    # This stuff actually matters.
    # Need the fluid velocity to be zero so we can ignore the denominator.
    y1dot = @SVector [0.0, 0.0, 0.0]
    y1dot_fluid = @SVector [0.0, 0.0, 0.0]
    y0dot = @SVector [0.0, 0.0, 0.0]
    chord_uvec = [1.0, 0.0, 0.0]
    span_uvec = [0.0, 1.0, 0.0]

    # Create a source element.
    se = AcousticAnalogies.TBLTESourceElement(c0, nu, dr, chord, y0dot, y1dot, y1dot_fluid, τ, dτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_suction_uvec)

    # Now, create an observer at different places, and check that we get the correct directivity function.
    for x_er in [-5.0, -3.0, 1.5, 4.0]
        for y_er in [-5.0, -3.0, 1.5, 4.0]
            for z_er in [-5.0, -3.0, 1.5, 4.0]
                x = @SVector [x_er, y_er, z_er]
                obs = AcousticAnalogies.StationaryAcousticObserver(x)
                r_er = sqrt(x_er^2 + y_er^2 + z_er^2)
                Θ_er = acos(x_er/r_er)
                Φ_er = acos(y_er/sqrt(y_er^2 + z_er^2))

                # Observer time doesn't matter since the observer is stationary.
                t_obs = 7.0
                Dl = (sin(Θ_er)^2) * (sin(Φ_er)^2)
                @test AcousticAnalogies.Dbar_l(se, obs, t_obs) ≈ Dl
                Dh = 2*(sin(0.5*Θ_er)^2) * (sin(Φ_er)^2)
                @test AcousticAnalogies.Dbar_h(se, obs, t_obs) ≈ Dh

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

                @test AcousticAnalogies.Dbar_l(se_trans, obs_trans, t_obs) ≈ Dl
                @test AcousticAnalogies.Dbar_h(se_trans, obs_trans, t_obs) ≈ Dh
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
            for vr in [0.0, 2.1, -2.1]
                for θ in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                    for twist in [5, 10, 65, 95, 260, 270, 290].*(pi/180)
                        for vn_sign in [1, -1]
                            for vc_sign in [1, -1]
                                vn = -(Vx + u)*vn_sign
                                vc = (-Vy + v)*vc_sign
                                se = AcousticAnalogies.TBLTESourceElement(c0, nu, r, θ, dr, chord, twist, vn, vr, vc, τ, dτ, bl, twist_about_positive_y)
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

end # module
