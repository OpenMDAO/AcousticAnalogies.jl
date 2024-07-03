module BPMITRTests

using AcousticAnalogies
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: Akima
using StaticArrays: @SVector
using Test

const aspect_data = [2.0,2.67,4.0,6.0,12.0,24.0]
const aratio_data = [0.54,0.62,0.71,0.79,0.89,0.95]
const aspect_ratio_correction = Akima(aspect_data, aratio_data)

struct BMTipAlphaCorrection{TCorrection} <: AcousticAnalogies.AbstractTipAlphaCorrection
    correction::TCorrection

    function BMTipAlphaCorrection(aspect_ratio)
        # Taken from BYUFLOWLab's BPM.jl package.
        if aspect_ratio < 2.0
            correction = 0.5*one(aspect_ratio)
        elseif aspect_ratio <= 24.0
            correction = aspect_ratio_correction(aspect_ratio)
        else
            correction = 1.0*one(aspect_ratio)
        end
        return new{typeof(correction)}(correction)
    end
end

function AcousticAnalogies.tip_vortex_alpha_correction(blade_tip::AcousticAnalogies.AbstractBladeTip{<:BMTipAlphaCorrection}, alphatip)
    return AcousticAnalogies.tip_alpha_correction(blade_tip).correction*alphatip
end

@testset "figure 22b" begin
    data_input = load(joinpath(@__DIR__, "gen_bpmjl_data", "figure22b-inputs.jld2"))
    rho = data_input["rho"]
    asound = data_input["asound"]
    mu = data_input["mu"]
    Vinf = data_input["Vinf"]
    omega = data_input["omega"]
    B = data_input["B"]
    Rhub = data_input["Rhub"]
    Rtip = data_input["Rtip"]
    radii = data_input["radii"]
    chord = data_input["chord"]
    twist = data_input["twist"]
    alpha = data_input["alpha"]
    U = data_input["U"]
    hs = data_input["hs"]
    Psis = data_input["Psis"]
    num_src_times_blade_pass = data_input["num_src_times_blade_pass"]
    tripped_flags = data_input["tripped_flags"]

    num_radial = length(radii)

    nu = mu/rho
    dradii = AcousticAnalogies.get_dradii(radii, Rhub, Rtip)

    # Get some transform stuff.
    bpp = 1/(B/(2*pi)*omega)  # 1/(B blade_passes/rev * 1 rev / (2*pi rad) * omega rad/s)
    num_blade_pass = 1
    period_src = num_blade_pass*bpp
    num_src_times = num_src_times_blade_pass * num_blade_pass
    t0 = 0.0
    dt = period_src/num_src_times
    src_times = t0 .+ (0:num_src_times-1).*dt

    # I don't see any discussion for what type of tip was used for the tip vortex noise.
    # FlatTip with no CCBlade.jl tip correction  or BPM-style tip correction seems to match the BARC predictions well.
    # blade_tip = AcousticAnalogies.FlatTip(AcousticAnalogies.NoTipAlphaCorrection())
    # BPM.jl uses a different tip alpha correction which appears to require the blade aspect ratio.
    # First need the blade area.
    area = sum(chord .* dradii)
    cbar = area/(Rtip - Rhub)
    aspect_ratio = Rtip/cbar
    alpha0lift = 0.0
    blade_tip = AcousticAnalogies.FlatTip(BMTipAlphaCorrection(aspect_ratio), alpha0lift)

    # Getting the coordinate system consistent with BPM.jl is a bit tricky.
    # Here's a bit of code from BPM.jl:
    #
    # # Calculate the trailing edge position relative to the hub
    # xs = sin(beta)*d - cos(beta)*(c - c1)
    # zs = cos(beta)*d + sin(beta)*(c - c1)
    #
    # OK, so that shows me that the blade is initially aligned with the z axis, rotating to the positive x direction.
    # And I know the blades are rotating about the positive y axis.
    # So that's the answer for the BPM.jl coordinate system:
    #
    #   * freestream in the positive y axis direction.
    #   * first blade initially aligned with the positive z axis, rotating about the positive y axis.
    #
    # Now, what do I need to do with AcousticAnalogies to make that happen?
    # I want the blades to be translating in the negative y direction, rotating about the positive y axis.
    # I usually start with the blades rotating about either the positive or negative x axis, moving in the direction of the positive x axis.
    # I think the answer is,
    #
    #   * start out with the blades rotating about the negative x axis, moving in the direction of the positive x axis
    #   * rotate 90° about the negative z axis.
    #     After this, the blades will be moving in the negative y direction, rotating about the positive y axis, which is good.
    #     But I want the first blade to be aligned with the positive z axis, and stopping here would mean it's aligned with the positive x axis.
    #   * rotate 90° about the negative y axis.
    #     This will put the first blade in line with the positive z axis.

    # So, let's do what we said we need to do.
    # Start with a rotation about the negative x axis.
    positive_x_rotation = false
    rot_trans = SteadyRotXTransformation(t0, omega*ifelse(positive_x_rotation, 1, -1), 0)

    # Then translate along the positive x axis.
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    v0_hub = @SVector [Vinf, 0.0, 0.0]
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # Then a 90° rotation about the negative z axis.
    trans_z90deg = SteadyRotZTransformation(0.0, 0.0, -0.5*pi)

    # Then a 90° rotation about the negative y axis.
    trans_y90deg = SteadyRotYTransformation(0.0, 0.0, -0.5*pi)

    # Put them all together:
    trans = compose.(src_times, Ref(trans_y90deg),
                compose.(src_times, Ref(trans_z90deg),
                    compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))))

    # Use the M_c = 0.8*M that BPM.jl and the BPM report use.
    U = @. 0.8*sqrt(Vinf^2 + (omega*radii)^2)

    θs = (0:(B-1)) .* (2*pi/B) .* ifelse(positive_x_rotation, 1, -1)

    bl = [ifelse(tf,
                 AcousticAnalogies.TrippedN0012BoundaryLayer(),
                 AcousticAnalogies.UntrippedN0012BoundaryLayer()) for tf in tripped_flags]
    @show bl

    # Do it.
    
    θs_rs = reshape(θs, 1, 1, :)
    radii_rs = reshape(radii, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    chord_rs = reshape(chord, 1, :, 1)
    twist_rs = reshape(twist, 1, :, 1)
    hs_rs = reshape(hs, 1, :, 1)
    Psis_rs = reshape(Psis, 1, :, 1)
    Us_rs = reshape(U, 1, :, 1)
    alphas_rs = reshape(alpha, 1, :, 1)
    bls_rs = reshape(bl, 1, :, 1)

    # Separate things into tip and no-tip.
    radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
    dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
    chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
    twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
    hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
    Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
    Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
    alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]
    bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

    radii_rs_with_tip = @view radii_rs[:, end:end, :]
    dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
    chord_rs_with_tip = @view chord_rs[:, end:end, :]
    twist_rs_with_tip = @view twist_rs[:, end:end, :]
    hs_rs_with_tip = @view hs_rs[:, end:end, :]
    Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
    Us_rs_with_tip = @view Us_rs[:, end:end, :]
    alphas_rs_with_tip = @view alphas_rs[:, end:end, :]
    bls_rs_with_tip = @view bls_rs[:, end:end, :]

    direct = AcousticAnalogies.BPMDirectivity
    use_UInduction = false
    use_Doppler = false
    mach_correction = AcousticAnalogies.NoMachCorrection
    ses_no_tip = CombinedNoTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord_rs_no_tip, twist_rs_no_tip, hs_rs_no_tip, Psis_rs_no_tip, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans

    ses_with_tip = CombinedWithTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord_rs_with_tip, twist_rs_with_tip, hs_rs_with_tip, Psis_rs_with_tip, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

    # Write out the source elements.
    pvd_no_tip = AcousticAnalogies.to_paraview_collection(joinpath(@__DIR__, "figure22b-no_tip"), ses_no_tip)
    pvd_with_tip = AcousticAnalogies.to_paraview_collection(joinpath(@__DIR__, "figure22b-with_tip"), ses_with_tip)

    # Put the source elements together:
    ses = cat(ses_no_tip, ses_with_tip; dims=2)

end


end # module
