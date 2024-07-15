using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: Akima
using StaticArrays
using Test

data = load(joinpath(@__DIR__, "gen_bpmjl_data", "figure22b.jld2"))
rho = data["rho"]
asound = data["asound"]
mu = data["mu"]
Vinf = data["Vinf"]
omega = data["omega"]
B = data["B"]
Rhub = data["Rhub"]
Rtip = data["Rtip"]
radii = data["radii"]
chord = data["chord"]
twist = data["twist"]
alpha = data["alpha"]
U = data["U"]
hs = data["hs"]
Psis = data["Psis"]
num_src_times_blade_pass = data["num_src_times_blade_pass"]
tripped_flags = data["tripped_flags"]

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
# BPM.jl uses a different tip alpha correction which appears to require the blade aspect ratio, defined as the blade radius divided by the average chord.
cbar = sum(chord .* dradii) / (Rtip - Rhub)
aspect_ratio = Rtip/cbar
alpha0lift = 0.0
blade_tip = AcousticAnalogies.FlatTip(AcousticAnalogies.BMTipAlphaCorrection(aspect_ratio), alpha0lift)

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

# Azimuthal offset for each blade.
θs = (0:(B-1)) .* (2*pi/B) .* ifelse(positive_x_rotation, 1, -1)

bls = [ifelse(tf,
              AcousticAnalogies.TrippedN0012BoundaryLayer(),
              AcousticAnalogies.UntrippedN0012BoundaryLayer()) for tf in tripped_flags]

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing.
# bls_lblvs = fill(AcousticAnalogies.UntrippedN0012BoundaryLayer(), num_radial)
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

r_obs = 2.27 # meters
theta_obs = -35*pi/180
# So, the docstring for BPM.jl says that `V` argument is the wind velocity in the y direction.
# So I guess we should assume that the blades are rotating about the y axis.
# And if the freestream velocity is in the positive y axis, then, from the perspective of the fluid, the blades are translating in the negative y direction.
# And I want the observer to be downstream/behind the blades, so that would mean they would have a positive y position.
# So I want to rotate the observer around the positive x axis, so I'm going to switch the sign on `theta_obs`.
t0_obs = 0.0
x0_obs = [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]
# The observer is moving in the same direction as the blades, which is the negative y axis.
v_obs = @SVector [0.0, -Vinf, 0.0]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0_obs, x0_obs, v_obs)

# Reshape the inputs to the source element constructors so that everything will line up with (num_times, num_radial, num_blades).
θs_rs = reshape(θs, 1, 1, :)
radii_rs = reshape(radii, 1, :, 1)
dradii_rs = reshape(dradii, 1, :, 1)
chord_rs = reshape(chord, 1, :, 1)
twist_rs = reshape(twist, 1, :, 1)
hs_rs = reshape(hs, 1, :, 1)
Psis_rs = reshape(Psis, 1, :, 1)
Us_rs = reshape(U, 1, :, 1)
alphas_rs = reshape(alpha, 1, :, 1)
bls_rs = reshape(bls, 1, :, 1)
# bls_untripped_rs = reshape(bls_lblvs, 1, :, 1)

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

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing.
lblvs_ses_untripped = AcousticAnalogies.LBLVSSourceElement{direct,use_UInduction,use_Doppler}.(asound, nu, radii_rs, θs_rs, dradii_rs, chord_rs, twist_rs, Us_rs, alphas_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# # Write out the source elements.
# pvd_no_tip = AcousticAnalogies.to_paraview_collection(joinpath(@__DIR__, "figure22b-no_tip"), ses_no_tip)
# pvd_with_tip = AcousticAnalogies.to_paraview_collection(joinpath(@__DIR__, "figure22b-with_tip"), ses_with_tip)
# pvd_all = AcousticAnalogies.to_paraview_collection(joinpath(@__DIR__, "figure22b-all"), (ses_no_tip, ses_with_tip); observers=(obs,))

# Put the source elements together:
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Define the frequencies we'd like to evaluate.
# BPM.jl uses the approximate 1/3rd-octave bands.
freqs_obs = AcousticMetrics.ApproximateThirdOctaveCenterBands(100.0, 40000.0)
freqs_src = freqs_obs

# Now do the noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))
pbs_lblvss_untripped = AcousticAnalogies.noise.(lblvs_ses_untripped, Ref(obs), Ref(freqs_src))

# Separate out each source.
pbs_tblte_ps = AcousticAnalogies.pbs_pressure.(bpm_outs)
pbs_tblte_ss = AcousticAnalogies.pbs_suction.(bpm_outs)
pbs_tblte_alphas = AcousticAnalogies.pbs_alpha.(bpm_outs)
pbs_tebs = AcousticAnalogies.pbs_teb.(bpm_outs)
pbs_tips = AcousticAnalogies.pbs_tip.(bpm_outs[:, end:end, :])

# Combine each noise prediction.
time_axis = 1
pbs_pressure = AcousticMetrics.combine(pbs_tblte_ps, freqs_obs, time_axis)
pbs_suction = AcousticMetrics.combine(pbs_tblte_ss, freqs_obs, time_axis)
pbs_alpha = AcousticMetrics.combine(pbs_tblte_alphas, freqs_obs, time_axis)
pbs_teb = AcousticMetrics.combine(pbs_tebs, freqs_obs, time_axis)
pbs_tip = AcousticMetrics.combine(pbs_tips, freqs_obs, time_axis)
pbs_lblvs_untripped = AcousticMetrics.combine(pbs_lblvss_untripped, freqs_obs, time_axis)

# Now I need to account for the fact that Figure 22b is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
# So, to do that, I need to multiply the mean-squared pressure by Δf_nb/Δf_pbs, where `Δf_nb` is the 20 Hz narrowband and `Δf_pbs` is the bandwidth of each 1/3-octave proportional band.
# (Dividing the MSP by Δf_pbs aka the 1/3 octave spacing is like getting a power-spectral density, then multiplying by the narrowband spacing Δf_nb gives us the MSP associated with the narrowband.)
# I think the paper describes that, right?
# Right, here's something:
#
#   > The current prediction method is limited to one-third octave bands, but it is compared to the narrowband experiment with Δf = 20 Hz.
#   > This is done by dividing the energy from the one-third octave bands by the number of bands in Δf = 20 Hz.
#
# So, `Δf_pbs/Δf_nb` would represent the number of `Δf_nb`-width bands that could fit in a proportional band of bin width `Δf_pbs`.
# And then I'm dividing by that.
# So that seems like the right thing.
# So, first thing is to get the proportional band spacing.
freqs_l = AcousticMetrics.lower_bands(freqs_obs)
freqs_u = AcousticMetrics.upper_bands(freqs_obs)
df_pbs = freqs_u .- freqs_l
# Also need the experimental narrowband spacing.
df_nb = 20.0
# Now multiply each by that.
nb_pressure = pbs_pressure .* df_nb ./ df_pbs
nb_suction = pbs_suction .* df_nb ./ df_pbs
nb_alpha = pbs_alpha .* df_nb ./ df_pbs
nb_teb = pbs_teb .* df_nb ./ df_pbs
nb_tip = pbs_tip .* df_nb ./ df_pbs
nb_lblvs_untripped = pbs_lblvs_untripped .* df_nb ./ df_pbs

# Now I want the SPL, which should just be this:
pref = 20e-6
spl_pressure = 10 .* log10.(nb_pressure./(pref^2))
spl_suction = 10 .* log10.(nb_suction./(pref^2))
spl_alpha = 10 .* log10.(nb_alpha./(pref^2))
spl_teb = 10 .* log10.(nb_teb./(pref^2))
spl_tip = 10 .* log10.(nb_tip./(pref^2))
spl_lblvs_untripped = 10 .* log10.(nb_lblvs_untripped./(pref^2))

# Read in the BPM.jl data.
freq_bpmjl = data["freqs"]
spl_pressure_bpmjl = data["spl_nb_pressure"]
spl_suction_bpmjl = data["spl_nb_suction"]
spl_separation_bpmjl = data["spl_nb_separation"]
spl_lblvs_bpmjl = data["spl_nb_lblvs"]
spl_blunt_bpmjl = data["spl_nb_blunt"]
spl_tip_bpmjl = data["spl_nb_tip"]

# The frequencies in the CSV file should match the observer frequencies we're using.
@test all(freqs_obs .≈ freq_bpmjl)

# Only look at the SPLs that are actually significant, i.e. greater than 10 dB.
@test maximum(abs.(spl_pressure[spl_pressure_bpmjl .> 10] .- spl_pressure_bpmjl[spl_pressure_bpmjl .> 10])) < 0.77
@test maximum(abs.(spl_suction[spl_suction_bpmjl .> 10] .- spl_suction_bpmjl[spl_suction_bpmjl .> 10])) < 0.78
@test maximum(abs.(spl_alpha[spl_separation_bpmjl .> 10] .- spl_separation_bpmjl[spl_separation_bpmjl .> 10])) < 0.78
@test maximum(abs.(spl_lblvs_untripped[spl_lblvs_bpmjl .> 10] .- spl_lblvs_bpmjl[spl_lblvs_bpmjl .> 10])) < 0.81
@test maximum(abs.(spl_teb[spl_blunt_bpmjl .> 10] .- spl_blunt_bpmjl[spl_blunt_bpmjl .> 10])) < 0.81
@test maximum(abs.(spl_tip[spl_tip_bpmjl .> 10] .- spl_tip_bpmjl[spl_tip_bpmjl .> 10])) < 1.25

