module AdvancedTimeTests

using AcousticAnalogies
using LinearAlgebra: norm
using NLsolve
using StaticArrays
using Test

@testset "Advanced time tests" begin
    # Goal is to verify that the code can solve the equation
    #
    #   R(t) = t - (τ + |x(t) - y|/c0) = 0
    #
    # for t, where τ and y are the source time and position, t and x are the
    # observer time and position, and c0 is the (constant) speed of sound.

    # Create a source element for the test.
    # The only things about the source element that matters to the advanced
    # time calculation is the time and position, and the speed of sound.
    τ = 2.5
    y = @SVector [-4.0, 3.0, 6.0]
    c0 = 2.0
    dummy0 = 1.0
    dummy3 = @SVector [0.0, 0.0, 0.0]
    se = CompactSourceElement(dummy0, c0, dummy0, dummy0, y, dummy3, dummy3, dummy3, dummy3, dummy3, τ, dummy3)

    # Define a function that takes a source element and an observer and compares
    # the advanced time solution to one found by NLsolve.
    function compare_to_nlsolve(se, obs)
        # Create the residual equation that we'll solve. nlsolve assumes the
        # residual function takes in and returns arrays.
        R(t) = [t[1] - (se.τ + norm(obs(t[1]) .- se.y0dot)/se.c0)]

        # Solve the advanced time equation.
        result = nlsolve(R, [1.0], autodiff=:forward)
        if !converged(result)
            @error "nlsolve advanced time calculation did not converge:\n$(result)"
        end
        t_obs = result.zero[1]

        # Check that we can get the right answer.
        @test adv_time(se, obs) ≈ t_obs

        return nothing
    end

    @testset "Stationary observer" begin
        x0 = @SVector [-3.0, 2.0, 8.5]
        obs = StationaryAcousticObserver(x0)
        compare_to_nlsolve(se, obs)
    end

    @testset "Constant velocity observer" begin
        t0 = 3.5
        x0 = @SVector [-2.0, 3.5, 6.25]
        v = @SVector [-1.5, 1.5, 3.25]
        obs = ConstVelocityAcousticObserver(t0, x0, v)
        compare_to_nlsolve(se, obs)
    end

end


end # module
