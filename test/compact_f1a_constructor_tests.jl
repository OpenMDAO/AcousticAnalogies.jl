module CompactF1AConstructorTests

using AcousticAnalogies
using CCBlade
using DelimitedFiles
using KinematicCoordinateTransformations
using StaticArrays
using JLD2
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
end

@testset "CCBlade CompactSourceElement complete test" begin

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
