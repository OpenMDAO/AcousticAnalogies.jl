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
    Vinf = ccbc.v
    omega = 2.0
    precone = 3*pi/180
    rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; precone=precone, turbine=false)
    ops = simple_op.(Vinf, omega, r, ccbc.rho; precone=precone)

    dummy = similar(r)
    sections = Section.(r, dummy, dummy, dummy)

    @test AcousticAnalogies.get_ccblade_omega(rotor, sections, ops) ≈ omega
    @test AcousticAnalogies.get_ccblade_Vinf(rotor, sections, ops) ≈ Vinf
    @test all(AcousticAnalogies.get_ccblade_dradii(rotor, sections) .≈ Δr)
end

@testset "CCBlade CompactSourceElement test" begin

    omega = 2200*(2*pi/60)
    # Read in the loading data.
    data = readdlm("gen_test_data/gen_ccblade_data/ccblade_omega11.csv", ',')
    fn = data[:, 1]
    fc = data[:, 2]

    # Create the CCBlade objects.
    rotor = Rotor(ccbc.Rhub, ccbc.Rtip, ccbc.num_blades; turbine=false)
    sections = Section.(ccbc.radii, ccbc.chord, ccbc.theta, nothing)
    ops = simple_op.(ccbc.v, omega, ccbc.radii, ccbc.rho; asound=ccbc.c0)
    # Only care about getting the loading into the CCBlade Output structs.
    dummies = fill(0.0, 13)
    outs = Outputs.(fn, fc, dummies...)

    # Finally get all the source elements.
    num_blade_passes = 3
    steps_per_blade_pass = 8
    ses_helper = source_elements_ccblade(rotor, sections, ops, outs, ccbc.area_over_chord_squared, num_blade_passes, steps_per_blade_pass)

    # Now need to get the source elements the "normal" way. First get the
    # transformation objects.
    rot_axis = @SVector [1.0, 0.0, 0.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    v0_hub = ccbc.v.*rot_axis
    t0 = 0.0
    rot_trans = SteadyRotXTransformation(t0, omega, 0.0)
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # Need the source times.
    # Blade passing period (amount of time for one blade pass).
    bpp = 2*pi/omega/ccbc.num_blades
    src_time_range = num_blade_passes*bpp
    num_src_times = num_blade_passes*steps_per_blade_pass
    dt = src_time_range/(num_src_times - 1)
    src_times = t0 .+ (0:num_src_times-1).*dt

    # This is just an array of the angular offsets of each blade.
    θs = 2*pi/ccbc.num_blades.*(0:(ccbc.num_blades-1))

    # Radial spacing.
    dradii = get_dradii(ccbc.radii, ccbc.Rhub, ccbc.Rtip)

    # Reshape stuff for broadcasting.
    radii = reshape(ccbc.radii, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    cs_area = reshape(ccbc.area_over_chord_squared.*ccbc.chord.^2, 1, :, 1)
    fn = reshape(fn, 1, :, 1)
    fc = reshape(fc, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This isn't really necessary.
    θs = reshape(θs, 1, 1, :)

    # Get all the transformations.
    trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

    # Transform the source elements.
    ses = CompactSourceElement.(ccbc.rho, ccbc.c0, radii, θs, dradii, cs_area, fn, fc, src_times) .|> trans

    for field in fieldnames(CompactSourceElement)
        @test all(getproperty.(ses_helper, field) .≈ getproperty.(ses, field))
    end
end

end
