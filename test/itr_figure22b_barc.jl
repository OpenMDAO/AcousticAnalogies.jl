using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: readdlm
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: akima
using StaticArrays: @SVector
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
num_src_times_blade_pass = data["num_src_times_blade_pass"]
Psis = data["Psis"]

nu = mu/rho
num_radial = length(radii)
dradii = AcousticAnalogies.get_dradii(radii, Rhub, Rtip)

# Get the source time, which will be one blade pass worth of time, each blade pass with `num_src_times_blade_pass` steps per blade pass.
bpp = 1/(B/(2*pi)*omega)  # 1/(B blade_passes/rev * 1 rev / (2*pi rad) * omega rad/s)
num_blade_pass = 1
period_src = num_blade_pass*bpp
num_src_times = num_src_times_blade_pass * num_blade_pass
t0 = 0.0
dt = period_src/num_src_times
src_times = t0 .+ (0:num_src_times-1).*dt

# Now let's define the coordinate system.
# I'm going to do my usual thing, which is to have the freestream velocity pointed in the negative x direction, and thus the blades will be translating in the positive x direction.
# And the blades will be rotating about the positive x axis at a rate of `omega`.
rot_trans = SteadyRotXTransformation(t0, omega, 0.0)

# The hub/rotation axis of the blades will start at the origin at time `t0`, and translate in the positive x direction at a speed of `Vinf`.
y0_hub = @SVector [0.0, 0.0, 0.0]  # m
v0_hub = @SVector [Vinf, 0.0, 0.0] # m/s
const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

# Now I can put the two transformations together:
trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

# Azimuthal offset for each blade.
θs = (0:(B-1)) .* (2*pi/B)

# Paper doesn't specify the microphone used for Figure 22, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 22.
# For the coordinate system, I'm doing my usual thing, which is to have the freestream velocity pointed in the negative x direction, and thus the blades will be translating in the positive x direction.
# The observer (microphone 5) is 35 deg behind/downstream of the rotor rotation plane, so this should be good.
# But it will of course be moving with the same freestream in the positive x direction.
r_obs = 2.27 # meters
theta_obs = -35*pi/180
# The observer is moving in the positive x direction at Vinf, at the origin at time t0.
t0_obs = 0.0
x0_obs = @SVector [r_obs*sin(theta_obs), r_obs*cos(theta_obs), 0.0]
v_obs = @SVector [Vinf, 0.0, 0.0]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0_obs, x0_obs, v_obs)

# In the text describing Figure 22, "For these predictions, the trip flag was set to “tripped”, due to the rough surface quality of the blade."
bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
# In the Figure 22 caption, "for these predictions, bluntness thickness H was set to 0.8 mm and trailing edge angle Ψ was set to 16 degrees."
h = 0.8e-3  # meters
Psi = 16*pi/180  # radians

# I don't see any discussion for what type of tip was used for the tip vortex noise.
# The flat tip seems to match the PAS+ROTONET+BARC predictions well.
blade_tip = AcousticAnalogies.FlatTip()

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
# bls_rs = reshape(bls, 1, :, 1)

# Separate things into tip and no-tip.
radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]
# bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

radii_rs_with_tip = @view radii_rs[:, end:end, :]
dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
chord_rs_with_tip = @view chord_rs[:, end:end, :]
twist_rs_with_tip = @view twist_rs[:, end:end, :]
hs_rs_with_tip = @view hs_rs[:, end:end, :]
Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
Us_rs_with_tip = @view Us_rs[:, end:end, :]
alphas_rs_with_tip = @view alphas_rs[:, end:end, :]
# bls_rs_with_tip = @view bls_rs[:, end:end, :]

positive_x_rotation = true
ses_no_tip = CombinedNoTipBroadbandSourceElement.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord_rs_no_tip, twist_rs_no_tip, hs_rs_no_tip, Psis_rs_no_tip, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, Ref(bl), positive_x_rotation) .|> trans

ses_with_tip = CombinedWithTipBroadbandSourceElement.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord_rs_with_tip, twist_rs_with_tip, hs_rs_with_tip, Psis_rs_with_tip, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, Ref(bl), Ref(blade_tip), positive_x_rotation) .|> trans

# It's more convinient to cat all the sources together.
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# The predictions in Figure 22b appear to be on 1/3 octave, ranging from about 200 Hz to 60,000 Hz.
# But let's expand the range of source frequencies to account for Doppler shifting.
freqs_src = AcousticMetrics.ExactProportionalBands{3, :center}(10.0, 200000.0)
freqs_obs = AcousticMetrics.ExactProportionalBands{3, :center}(200.0, 60000.0)

# Now we can do a noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))

# This seperates out the noise prediction for each source-observer combination into the different sources.
pbs_tblte_ps = AcousticAnalogies.pbs_pressure.(bpm_outs)
pbs_tblte_ss = AcousticAnalogies.pbs_suction.(bpm_outs)
pbs_tblte_alphas = AcousticAnalogies.pbs_alpha.(bpm_outs)
pbs_lblvss = AcousticAnalogies.pbs_lblvs.(bpm_outs)
pbs_tebs = AcousticAnalogies.pbs_teb.(bpm_outs)
pbs_tips = AcousticAnalogies.pbs_tip.(bpm_outs[:, end:end, :])

# Now, need to combine each broadband noise prediction.
# The time axis the axis over which the time varies for each source.
time_axis = 1
pbs_pressure = AcousticMetrics.combine(pbs_tblte_ps, freqs_obs, time_axis)
pbs_suction = AcousticMetrics.combine(pbs_tblte_ss, freqs_obs, time_axis)
pbs_alpha = AcousticMetrics.combine(pbs_tblte_alphas, freqs_obs, time_axis)
pbs_lblvs = AcousticMetrics.combine(pbs_lblvss, freqs_obs, time_axis)
pbs_teb = AcousticMetrics.combine(pbs_tebs, freqs_obs, time_axis)
pbs_tip = AcousticMetrics.combine(pbs_tips, freqs_obs, time_axis)

# Now I need to account for the fact that Figure 22b is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
# So, to do that, I need to multiply the mean-squared pressure by Δf_nb/Δf_pbs, where `Δf_nb` is the 20 Hz narrowband and `Δf_pbs` is the bandwidth of each 1/3-octave proportional band.
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
nb_lblvs = pbs_lblvs .* df_nb ./ df_pbs
nb_teb = pbs_teb .* df_nb ./ df_pbs
nb_tip = pbs_tip .* df_nb ./ df_pbs

# Now I want the SPL, which should just be this:
pref = 20e-6
spl_pressure = 10 .* log10.(nb_pressure./(pref^2))
spl_suction = 10 .* log10.(nb_suction./(pref^2))
spl_alpha = 10 .* log10.(nb_alpha./(pref^2))
spl_lblvs = 10 .* log10.(nb_lblvs./(pref^2))
spl_teb = 10 .* log10.(nb_teb./(pref^2))
spl_tip = 10 .* log10.(nb_tip./(pref^2))

# Now I should be able to compare to the BARC data.
# Need to read it in first.
data_pressure_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-TBL-TE-pressure-2.csv"), ',')
freq_pressure_barc = data_pressure_barc[:, 1]
spl_pressure_barc = data_pressure_barc[:, 2]

data_suction_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-TBL-TE-suction-2.csv"), ',')
freq_suction_barc = data_suction_barc[:, 1]
spl_suction_barc = data_suction_barc[:, 2]

data_separation_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-separation-2.csv"), ',')
freq_separation_barc = data_separation_barc[:, 1]
spl_separation_barc = data_separation_barc[:, 2]

data_lblvs_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-LBLVS.csv"), ',')
freq_lblvs_barc = data_lblvs_barc[:, 1]
spl_lblvs_barc = data_lblvs_barc[:, 2]

data_teb_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-BVS.csv"), ',')
freq_teb_barc = data_teb_barc[:, 1]
spl_teb_barc = data_teb_barc[:, 2]

data_tip_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-tip_vortex_shedding.csv"), ',')
freq_tip_barc = data_tip_barc[:, 1]
spl_tip_barc = data_tip_barc[:, 2]

# Interpolate the AcousticAnalogies.jl data onto the frequencies from the BARC CSV file.
spl_pressure_interp = akima(freqs_obs, spl_pressure, freq_pressure_barc)
spl_suction_interp = akima(freqs_obs, spl_suction, freq_suction_barc)
spl_separation_interp = akima(freqs_obs, spl_alpha, freq_separation_barc)
spl_lblvs_interp = akima(freqs_obs, spl_lblvs, freq_lblvs_barc)
spl_teb_interp = akima(freqs_obs, spl_teb, freq_teb_barc)
spl_tip_interp = akima(freqs_obs, spl_tip, freq_tip_barc)

# Now compare.
@test maximum(abs.(spl_pressure_interp .- spl_pressure_barc)) < 0.543

# Lower frequencies don't line up as well as higher.
# Not sure why.
@test abs(spl_suction_interp[1] - spl_suction_barc[1]) < 3.59
@test abs(spl_suction_interp[2] - spl_suction_barc[2]) < 2.71
@test abs(spl_suction_interp[3] - spl_suction_barc[3]) < 1.24
@test maximum(abs.(spl_suction_interp[4:end] .- spl_suction_barc[4:end])) < 0.462

# Lower frequencies don't line up as well as higher.
# Not sure why.
@test all(abs.(spl_separation_interp[1:12] .- spl_separation_barc[1:12]) .< [8.54, 7.63, 7.48, 6.82, 6.51, 6.60, 5.99, 4.87, 4.21, 2.58, 1.29, 0.299])
@test maximum(abs.(spl_separation_interp[13:end] .- spl_separation_barc[13:end])) < 0.466

@test all(abs.(spl_teb_interp .- spl_teb_barc) .< [0.134, 0.349, 0.424, 0.133, 0.187, 0.603, 0.244, 2.60])
@test all(abs.(spl_tip_interp .- spl_tip_barc) .< [0.0641, 0.0426, 0.170, 0.0923, 0.116, 0.231, 0.228, 0.160, 0.121])

