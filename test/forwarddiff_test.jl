module ForwardDiffTests

using AcousticAnalogies
using ForwardDiff
using KinematicCoordinateTransformations
using StaticArrays
using CCBlade
using AcousticMetrics
using Test
using LinearAlgebra


function guided_example(x)
    Rx::eltype(x) = 0.0
    Rhub = x[1]
    Rtip = x[2]
    
    num_blades = 2

    radii = [
        0.92904E-01, 0.11751,  0.15631,  0.20097,
        0.24792    , 0.29563,  0.34336,  0.39068,
        0.43727    , 0.48291,  0.52741,  0.57060,
        0.61234    , 0.65249,  0.69092,  0.72752,
        0.76218    , 0.79479,  0.82527,  0.85352,
        0.87947    , 0.90303,  0.92415,  0.94275,
        0.95880    , 0.97224,  0.98304,  0.99117,
        0.99660    , 0.99932] .* Rtip

    θs = 2pi / num_blades .* (0:(num_blades-1)) .+ Rx

    dradii = get_dradii(radii, Rhub, Rtip)

    rho = 1.226  # kg/m^3
    c0 = 340.0  # m/s

    v = x[3]
    omega = x[4]  # rad/s

    cs_area_over_chord_squared = 0.064
    chord = x[5:34] .* Rtip
    cs_area = cs_area_over_chord_squared .* chord.^2

    fn = [32.87810395677037, 99.05130471878633,  190.1751697055377,
         275.9492967565419,  358.14423433748146, 439.64679797145624,
         520.1002808148281,  599.1445046901513,  676.2358818769462,
         751.3409657831587,  824.2087672338118,  894.4465814696498,
         961.9015451678036,  1026.0112737521583, 1086.2610633094212,
         1141.4900032393818, 1190.3376703335655, 1230.8999662260915,
         1260.375390697363,  1275.354422403355,  1271.8827617273287,
         1245.9059108698596, 1193.9967137923225, 1113.9397490286995,
         1005.273267675585,  869.4101036003673,  709.8100230053759,
         532.1946243370355,  346.53986082379265, 180.66763939805125] .+ Rx

    fc = [26.09881302938423, 55.5216259955307, 75.84767780212506, 
        84.84509232798283, 89.73045068624886, 93.02999477395113, 
        95.4384273852926, 97.31647535460424, 98.81063179767507, 
        100.07617771995163, 101.17251941705561, 102.11543878532882,
        102.94453631586998, 103.63835661864168, 104.18877957193807, 
        104.51732850056433, 104.54735678589765, 104.1688287897138, 
        103.20319203286938, 101.46246817378582, 99.11692436681635, 
        96.49009546562475, 93.45834266417528, 89.49783586366624,
        83.87176811707455, 75.83190739325453, 64.88004605331857, 
        50.98243352318318, 34.85525518071079, 19.358679206883078] .+ Rx

    period = 2pi / omega
    num_src_times = 64
    dt = 2 * period / (num_src_times-1)
    src_times = (0:num_src_times-1) .* dt

    θs = reshape(θs, 1, 1, :)
    radii = reshape(radii, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    cs_area = reshape(cs_area, 1, :, 1)
    fn = reshape(fn, 1, :, 1)
    fc = reshape(fc, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This isn't really necessary.

    ses = CompactF1ASourceElement.(rho, c0, radii, θs, dradii, cs_area, -fn, 0.0, fc, src_times)

    t0 = 0.0  # Time at which the angle between the source and target coordinate systems is equal to offest.
    offset = 0.0  # Angular offset between the source and target cooridante systems at t0.
    rot_trans = SteadyRotXTransformation(t0, omega, offset)

    rot_axis = @SVector [0.0, 0.0, 1.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    global_trans = ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))

    y0_hub = @SVector [0.0, 0.0, 0.0]  # Position of the hub at time t0
    v0_hub = SVector{3}(v.*rot_axis)   # Constant velocity of the hub in the global reference frame
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    trans = compose.(src_times, Ref(const_vel_trans), compose.(src_times, Ref(global_trans), Ref(rot_trans)))

    ses = ses .|> trans

    x0 = @SVector [100*12*0.0254, 0.0, 0.0]  # 100 ft in meters
    obs = StationaryAcousticObserver(x0)

    obs_time = adv_time.(ses, Ref(obs))

    apth = noise.(ses, Ref(obs), obs_time)

    bpp = period/num_blades  # blade passing period
    obs_time_range = 2 * bpp
    num_obs_times = 128
    apth_total = combine(apth, obs_time_range, num_obs_times, 1)

    oaspl_from_apth = AcousticMetrics.OASPL(apth_total)
    nbs = AcousticMetrics.MSPSpectrumAmplitude(apth_total)
    oaspl_from_nbs = AcousticMetrics.OASPL(nbs)

    return vcat(oaspl_from_apth, oaspl_from_nbs)
end

### Just run the guided example through ForwardDiff without errors.
# May want to make more comprehensive later. For now pretty much everything 
# is affected by Dual numbers here, except for the loads are still Floats 
# for now.

@testset "ForwardDiff test" begin

    Rhub = 0.10
    Rtip = 1.1684
    v = 0.0  # m/s
    omega = 2200 * 2pi/60  # rad/s

    chord_normalized = [0.35044     , 0.28260     , 0.22105     , 0.17787     , 0.14760,
                        0.12567     , 0.10927     , 0.96661E-01 , 0.86742E-01 , 0.78783E-01 , 
                        0.72287E-01 , 0.66906E-01 , 0.62387E-01 , 0.58541E-01 , 0.55217E-01 , 
                        0.52290E-01 , 0.49645E-01 , 0.47176E-01 , 0.44772E-01 , 0.42326E-01 , 
                        0.39732E-01 , 0.36898E-01 , 0.33752E-01 , 0.30255E-01 , 0.26401E-01 ,
                        0.22217E-01 , 0.17765E-01 , 0.13147E-01 , 0.85683E-02 , 0.47397E-02]

    x0 = vcat(Rhub, Rtip, v, omega, chord_normalized)

    J = ForwardDiff.jacobian(guided_example, x0)
end

end #module
