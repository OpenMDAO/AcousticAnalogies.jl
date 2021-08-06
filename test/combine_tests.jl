module CombineTests

using AcousticAnalogies
using FLOWMath: akima, linear
using Random
using Test

@testset "Combine F1AOutput tests" begin

    # Goal is to verify that the code can faithfully combine two acoustic
    # pressures on different time "grids" onto a single common grid.
    fa(t) = sin(2*pi*t) + 0.2*cos(4*pi*(t-0.1))
    fb(t) = cos(6*pi*t) + 0.3*sin(8*pi*(t-0.2))

    n = 101
    t1 = collect(range(0.0, 1.0, length=n))
    dt = t1[2] - t1[1]
    # Add a bit of random noise to the time grid. Make sure that the amount of
    # randomness isn't large enough to make the time values non-monotonically
    # increasing (i.e., they don't overlap).
    t1 .+= 0.49.*dt.*(1 - 2*rand(size(t1)))

    t2 = collect(range(0.1, 1.1, length=n))
    dt = t2[2] - t2[1]
    t2 .+= 0.49.*dt.*(1 - 2*rand(size(t2)))

    # Now I need a bunch of acoustic pressures.
    apth1 = @. F1AOutput(t1, fa(t1), 2*fa(t1))
    apth2 = @. F1AOutput(t2, fb(t2), 3*fb(t2))

    # Calculate the "exact" answer by coming up with a common time, then
    # evaluating the test functions directly on the common time grid.
    period = 0.5
    n_out = 51
    t_start = max(t1[1], t2[1])
    t_common = t_start .+ (0:n_out-1).*(period/n_out)

    p_m = @. fa(t_common)+fb(t_common)
    p_d = @. 2*fa(t_common)+3*fb(t_common)

    apth_test = F1AAcousticPressure(p_m, p_d, step(t_common), first(t_common))

    function combine_test_axis1(f_interp)
        # Put all the acoustic pressures in one array.
        apth = hcat(apth1, apth2)

        # Combine.
        axis = 1
        apth_out = combine(apth, period, n_out, axis, f_interp=f_interp)

        # Now find the scaled absolute error between the two approaches.
        p_m_min, p_m_max = extrema(apth_test.p_m)
        err = @. abs(apth_out.p_m - apth_test.p_m)/(p_m_max - p_m_min)
        err_max = maximum(err)
        @test err_max < 0.01

        p_d_min, p_d_max = extrema(apth_test.p_d)
        err = @. abs(apth_out.p_d - apth_test.p_d)/(p_d_max - p_d_min)
        err_max = maximum(err)
        @test err_max < 0.01

        return nothing
    end

    function combine_test_axis2(f_interp)
        # Put all the acoustic pressures in one array.
        apth = permutedims(hcat(apth1, apth2))

        # Combine.
        axis = 2
        apth_out = combine(apth, period, n_out, axis, f_interp=f_interp)

        # Now find the scaled absolute error between the two approaches.
        p_m_min, p_m_max = extrema(apth_test.p_m)
        err = @. abs(apth_out.p_m - apth_test.p_m)/(p_m_max - p_m_min)
        err_max = maximum(err)
        @test err_max < 0.01

        p_d_min, p_d_max = extrema(apth_test.p_d)
        err = @. abs(apth_out.p_d - apth_test.p_d)/(p_d_max - p_d_min)
        err_max = maximum(err)
        @test err_max < 0.01

        return nothing
    end

    @testset "linear interpolation" begin
        combine_test_axis1(linear)
        combine_test_axis2(linear)
    end

    @testset "Akima spline interpolation" begin
        combine_test_axis1(akima)
        combine_test_axis2(akima)
    end
end

end # module
