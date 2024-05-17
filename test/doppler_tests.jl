module DopplerTests

using AcousticAnalogies: AcousticAnalogies
using KinematicCoordinateTransformations: KinematicTransformation, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation, compose
using StaticArrays: SVector
using Test

# struct DummyElement{Ty0dot,Ty1dot,Ttime,TSoS} <: AcousticAnalogies.AbstractCompactSourceElement
struct DummyElement{Ty0dot,Ty1dot,Ttime,TSoS} <: AcousticAnalogies.AbstractBroadbandSourceElement{AcousticAnalogies.BPMDirectivity,false, AcousticAnalogies.NoMachCorrection}
    # Source position and its time derivatives.
    y0dot::Ty0dot
    y1dot::Ty1dot
    # Source time.
    τ::Ttime
    # speed of sound
    c0::TSoS
end

"""
    (trans::KinematicTransformation)(se::DummyElement)

Transform the position and orientation of a source element according to the coordinate system transformation `trans`.
"""
function (trans::KinematicTransformation)(se::DummyElement)
    linear_only = false
    y0dot, y1dot = trans(se.τ, se.y0dot, se.y1dot, linear_only)

    return DummyElement(y0dot, y1dot, se.τ, se.c0)
end

@testset "Doppler shift" begin
    @testset "no motion" begin
        @testset "stationary observer" begin
            obs = AcousticAnalogies.StationaryAcousticObserver(SVector(1.2, 2.3, 3.4))
            se = DummyElement(SVector(3.0, 4.3, 5.0), SVector(0.0, 0.0, 0.0), 10.0, 340.0)
            @test AcousticAnalogies.doppler_factor(se, obs) ≈ 1.0
        end

        @testset "\"moving\" observer" begin
            # Should get the same thing with a `ConstVelocityAcousticObserver` that's not moving.
            obs = AcousticAnalogies.ConstVelocityAcousticObserver(8.0, SVector(1.2, 2.3, 3.4), SVector(0.0, 0.0, 0.0))
            se = DummyElement(SVector(3.0, 4.3, 5.0), SVector(0.0, 0.0, 0.0), 10.0, 340.0)
            @test AcousticAnalogies.doppler_factor(se, obs) ≈ 1.0
        end
    end

    @testset "moving source" begin
        @testset "stationary observer" begin
            for mach_vector in -0.9:0.1:0.9
                c0 = 343.0

                # Observer time doesn't matter since the observer isn't moving.
                t_obs = 6.5
                x_obs = SVector(1.2, 2.3, 3.4)
                obs = AcousticAnalogies.StationaryAcousticObserver(x_obs)

                y0dot = SVector(1.2, 2.3, -3.4)
                y1dot = SVector(0.0, 0.0, mach_vector*c0)
                τ = 10.0
                se = DummyElement(y0dot, y1dot, τ, c0)

                doppler_factor_expected = 1/(1 - mach_vector)
                @test AcousticAnalogies.doppler_factor(se, obs) ≈ doppler_factor_expected

                # Now, rotate and translate both the source and the observer.
                # The Doppler shift factor should be the same, assuming we don't change the motion of the source and observer.
                # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
                trans1 = SteadyRotXTransformation(t_obs, 0.0, 3.0*pi/180)
                trans2 = SteadyRotYTransformation(t_obs, 0.0, 4.0*pi/180)
                trans3 = SteadyRotZTransformation(t_obs, 0.0, 5.0*pi/180)
                x_trans = SVector(2.0, 3.0, 4.0)
                v_trans = SVector(0.0, 0.0, 0.0)
                # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
                trans4 = ConstantVelocityTransformation(t_obs, x_trans, v_trans)

                # Transform the source and observer.
                trans = compose(t_obs, trans4, compose(t_obs, trans3, compose(t_obs, trans2, trans1)))
                se_trans = trans(se)
                obs_trans = AcousticAnalogies.StationaryAcousticObserver(trans(t_obs, obs(t_obs)))

                # Now we should still get the same Doppler factor.
                @test AcousticAnalogies.doppler_factor(se_trans, obs_trans) ≈ doppler_factor_expected
            end
        end

        @testset "\"moving\" observer" begin
            for mach_vector in -0.9:0.1:0.9
                c0 = 343.0

                # Observer time doesn't matter since the observer isn't moving.
                t_obs = 6.5
                x_obs = SVector(1.2, 2.3, 3.4)
                v_obs = SVector(0.0, 0.0, 0.0)
                obs = AcousticAnalogies.ConstVelocityAcousticObserver(t_obs, x_obs, v_obs)

                y0dot = SVector(1.2, 2.3, -3.4)
                y1dot = SVector(0.0, 0.0, mach_vector*c0)
                τ = 10.0
                se = DummyElement(y0dot, y1dot, τ, c0)

                doppler_factor_expected = 1/(1 - mach_vector)
                @test AcousticAnalogies.doppler_factor(se, obs) ≈ doppler_factor_expected

                # Now, rotate and translate both the source and the observer.
                # The Doppler shift factor should be the same, assuming we don't change the motion of the source and observer.
                # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
                trans1 = SteadyRotXTransformation(t_obs, 0.0, 3.0*pi/180)
                trans2 = SteadyRotYTransformation(t_obs, 0.0, 4.0*pi/180)
                trans3 = SteadyRotZTransformation(t_obs, 0.0, 5.0*pi/180)
                x_trans = SVector(2.0, 3.0, 4.0)
                v_trans = SVector(0.0, 0.0, 0.0)
                # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
                trans4 = ConstantVelocityTransformation(t_obs, x_trans, v_trans)

                # Transform the source and observer.
                trans = compose(t_obs, trans4, compose(t_obs, trans3, compose(t_obs, trans2, trans1)))
                se_trans = trans(se)
                obs_trans = AcousticAnalogies.ConstVelocityAcousticObserver(t_obs, trans(t_obs, obs(t_obs), AcousticAnalogies.velocity(t_obs, obs))...)

                # Now we should still get the same Doppler factor.
                @test AcousticAnalogies.doppler_factor(se_trans, obs_trans) ≈ doppler_factor_expected
            end
        end

        @testset "actual moving observer" begin
            for mach_vector_obs in -0.9:0.1:0.9
                for mach_vector_src in -0.9:0.1:0.9
                    c0 = 343.0

                    t_obs = 10.0
                    # Need to move the observer farther from the source for cases where they're moving towards each other, since the Doppler factor will change if they cross and move past each other.
                    x_obs = SVector(1.2, 2.3, 340.0)
                    v_obs = SVector(0.0, 0.0, mach_vector_obs*c0)
                    obs = AcousticAnalogies.ConstVelocityAcousticObserver(t_obs, x_obs, v_obs)

                    y0dot = SVector(1.2, 2.3, -3.4)
                    y1dot = SVector(0.0, 0.0, mach_vector_src*c0)
                    τ = 10.0
                    se = DummyElement(y0dot, y1dot, τ, c0)

                    doppler_factor_expected = (1 - mach_vector_obs)/(1 - mach_vector_src)
                    @test AcousticAnalogies.doppler_factor(se, obs) ≈ doppler_factor_expected

                    # Now, rotate and translate both the source and the observer.
                    # The Doppler shift factor should be the same, assuming we don't change the motion of the source and observer.
                    # Time parameter for the steady rotations doesn't matter because the rotation rate is zero.
                    trans1 = SteadyRotXTransformation(t_obs, 0.0, 3.0*pi/180)
                    trans2 = SteadyRotYTransformation(t_obs, 0.0, 4.0*pi/180)
                    trans3 = SteadyRotZTransformation(t_obs, 0.0, 5.0*pi/180)
                    x_trans = SVector(2.0, 3.0, 4.0)
                    v_trans = SVector(0.0, 0.0, 0.0)
                    # Time parameter for the constant velocity transformations doesn't matter because the velocity is zero.
                    trans4 = ConstantVelocityTransformation(t_obs, x_trans, v_trans)

                    # Transform the source and observer.
                    trans = compose(t_obs, trans4, compose(t_obs, trans3, compose(t_obs, trans2, trans1)))
                    se_trans = trans(se)
                    obs_trans = AcousticAnalogies.ConstVelocityAcousticObserver(t_obs, trans(t_obs, obs(t_obs), AcousticAnalogies.velocity(t_obs, obs))...)

                    # Now we should still get the same Doppler factor.
                    @test AcousticAnalogies.doppler_factor(se_trans, obs_trans) ≈ doppler_factor_expected
                end
            end
        end
    end
end

end # module
