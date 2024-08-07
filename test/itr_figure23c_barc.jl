using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: readdlm
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: akima
using StaticArrays: @SVector
using Test

data = load(joinpath(@__DIR__, "gen_bpmjl_data", "figure23c.jld2"))
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

# Paper doesn't specify the microphone used for Figure 23, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 23.
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

# So, for the boundary layer, we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
num_untripped = Int(round(0.95*num_radial))
num_tripped = num_radial - num_untripped
bls_untripped = fill(AcousticAnalogies.UntrippedN0012BoundaryLayer(), num_untripped)
bls_tripped = fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), num_tripped)
bls = vcat(bls_untripped, bls_tripped)

# Now, the other trick: need to only include LBLVS noise for elements where the Reynolds number is < 160000.
# So, we need the Reynolds number for each section.
Re_c = @. U * chord / nu
# So now we want to extract the radial stations that meet that < 160000 condition.
low_Re_c = 160000
mask_low_Re_c = Re_c .< low_Re_c

# And we're also going to use the untripped boundary layer for the LBLVS source.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# In the Figure 23 caption, "for these predictions, bluntness thickness H was set to 0.5 mm and trailing edge angle Ψ was set to 14 degrees."
h = 0.5e-3  # meters
Psi = 14*pi/180  # radians

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
bls_rs = reshape(bls, 1, :, 1)

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

positive_x_rotation = true
ses_no_tip = CombinedNoTipBroadbandSourceElement.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord_rs_no_tip, twist_rs_no_tip, hs_rs_no_tip, Psis_rs_no_tip, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans
ses_with_tip = CombinedWithTipBroadbandSourceElement.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord_rs_with_tip, twist_rs_with_tip, hs_rs_with_tip, Psis_rs_with_tip, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

# It's more convinient to cat all the sources together.
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS stuff separately.
# Grab the parts of the inputs that correspond to the low Reynolds number stations.
radii_lblvs = @view radii[mask_low_Re_c]
dradii_lblvs = @view dradii[mask_low_Re_c]
chord_lblvs = @view chord[mask_low_Re_c]
twist_lblvs = @view twist[mask_low_Re_c]
hs_lblvs = @view hs[mask_low_Re_c]
Psis_lblvs = @view Psis[mask_low_Re_c]
Us_lblvs = @view U[mask_low_Re_c]
alphas_lblvs = @view alpha[mask_low_Re_c]

# And do the reshaping.
radii_lblvs_rs = reshape(radii_lblvs, 1, :, 1)
dradii_lblvs_rs = reshape(dradii_lblvs, 1, :, 1)
chord_lblvs_rs = reshape(chord_lblvs, 1, :, 1)
twist_lblvs_rs = reshape(twist_lblvs, 1, :, 1)
hs_lblvs_rs = reshape(hs_lblvs, 1, :, 1)
Psis_lblvs_rs = reshape(Psis_lblvs, 1, :, 1)
Us_lblvs_rs = reshape(Us_lblvs, 1, :, 1)
alphas_lblvs_rs = reshape(alphas_lblvs, 1, :, 1)

# Now we can create the source elements.
ses_lblvs = LBLVSSourceElement.(asound, nu, radii_lblvs_rs, θs_rs, dradii_lblvs_rs, chord_lblvs_rs, twist_lblvs_rs, Us_lblvs_rs, alphas_lblvs_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# The predictions in Figure 23c appear to be on 1/3 octave, ranging from about 200 Hz to 60,000 Hz.
# But let's expand the range of source frequencies to account for Doppler shifting.
freqs_src = AcousticMetrics.ExactProportionalBands{3, :center}(10.0, 200000.0)
freqs_obs = AcousticMetrics.ExactProportionalBands{3, :center}(200.0, 60000.0)

# Now we can do a noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))
pbs_lblvss = AcousticAnalogies.noise.(ses_lblvs, Ref(obs), Ref(freqs_src))

# This seperates out the noise prediction for each source-observer combination into the different sources.
pbs_tblte_ps = AcousticAnalogies.pbs_pressure.(bpm_outs)
pbs_tblte_ss = AcousticAnalogies.pbs_suction.(bpm_outs)
pbs_tblte_alphas = AcousticAnalogies.pbs_alpha.(bpm_outs)
pbs_tebs = AcousticAnalogies.pbs_teb.(bpm_outs)
pbs_tips = AcousticAnalogies.pbs_tip.(bpm_outs[:, end:end, :])

# Now, need to combine each broadband noise prediction.
# The time axis the axis over which the time varies for each source.
time_axis = 1
pbs_pressure = AcousticMetrics.combine(pbs_tblte_ps, freqs_obs, time_axis)
pbs_suction = AcousticMetrics.combine(pbs_tblte_ss, freqs_obs, time_axis)
pbs_alpha = AcousticMetrics.combine(pbs_tblte_alphas, freqs_obs, time_axis)
pbs_teb = AcousticMetrics.combine(pbs_tebs, freqs_obs, time_axis)
pbs_tip = AcousticMetrics.combine(pbs_tips, freqs_obs, time_axis)
pbs_lblvs = AcousticMetrics.combine(pbs_lblvss, freqs_obs, time_axis)

# Now I need to account for the fact that Figure 23c is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
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
data_pressure_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-TBL-TE-pressure.csv"), ',')
freq_pressure_barc = data_pressure_barc[:, 1]
spl_pressure_barc = data_pressure_barc[:, 2]

data_suction_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-TBL-TE-suction.csv"), ',')
freq_suction_barc = data_suction_barc[:, 1]
spl_suction_barc = data_suction_barc[:, 2]

data_separation_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-separation.csv"), ',')
freq_separation_barc = data_separation_barc[:, 1]
spl_separation_barc = data_separation_barc[:, 2]

data_lblvs_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-LBLVS.csv"), ',')
freq_lblvs_barc = data_lblvs_barc[:, 1]
spl_lblvs_barc = data_lblvs_barc[:, 2]

data_teb_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-BVS.csv"), ',')
freq_teb_barc = data_teb_barc[:, 1]
spl_teb_barc = data_teb_barc[:, 2]

data_tip_barc = readdlm(joinpath(@__DIR__, "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-tip_vortex_shedding.csv"), ',')
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
@test all(abs.(spl_pressure_interp .- spl_pressure_barc) .< [1.918, 1.880, 1.484, 1.65, 1.496, 1.170, 1.043, 0.729, 0.406, 0.406, 1.49, 1.142, 1.131])
@test all(abs.(spl_suction_interp .- spl_suction_barc) .< [2.193, 2.066, 1.984, 1.961, 1.686, 1.423, 1.255, 1.060, 0.339, 0.101, 0.149, 0.749, 1.363, 1.220, 1.547, 1.979])
@test all(abs.(spl_separation_interp .- spl_separation_barc) .< [17.002, 14.84, 12.09, 10.20, 9.42, 8.371, 7.763, 7.504, 7.099, 6.124, 5.307, 2.843, 2.326, 2.560, 2.583, 2.088, 1.448, 0.628, 0.112, 0.873, 1.971])
@test all(abs.(spl_lblvs_interp .- spl_lblvs_barc) .< [3.369, 3.795, 3.758, 3.797, 3.765, 3.749, 3.545, 3.927, 3.922, 3.652, 3.571])
@test all(abs.(spl_teb_interp .- spl_teb_barc) .< [0.274, 0.135, 0.211, 0.127, 0.0584, 2.200, 2.981])
@test all(abs.(spl_tip_interp .- spl_tip_barc) .< [0.590, 0.659, 0.625, 0.460, 0.240, 0.467, 0.434, 0.0235, 0.0468])
