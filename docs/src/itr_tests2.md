```@meta
CurrentModule = AADocs
```
# Software Quality Assurance, Cont.

## PAS/ROTONET/BARC Comparisons for the Pettingill et al. Ideally Twisted Rotor 

### Figure 22b
```@example figure22b
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: readdlm
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, ConstantVelocityTransformation
using FileIO: load
using GLMakie
using StaticArrays: @SVector

# Pettingill et al., "Acoustic And Performance Characteristics of an Ideally Twisted Rotor in Hover", 2021
# Parameters from Table 1
B = 4  # number of blades
Rtip = 0.1588 # meters
chord = 0.2*Rtip

# Standard day:
Tamb = 15 + 273.15 # 15°C in Kelvin
pamb = 101325.0  # Pa
R = 287.052874 # J/(kg*K)
rho = pamb/(R*Tamb)
asound = sqrt(1.4*R*Tamb)
# Dynamic and kinematic viscosity
mu = rho*1.4502e-5
nu = mu/rho

# This is a hover case, so the freestream velocity should be zero.
# CCBlade.jl will run with a zero freestream, but I've found that it compares a bit better with experiment if I give it a small non-zero value.
Vinf = 0.001*asound

# Figure 22 caption says Ω_c = 5465 RPM.
rpm = 5465.0
omega = rpm * (2*pi/60)

# Get "cell-centered" radial locations, and also the radial spacing.
num_radial = 50
r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
radii = r_Rtip .* Rtip
dradii = (r_Rtip_[2:end] .- r_Rtip_[1:end-1]) .* Rtip
Rhub = r_Rtip_[1]*Rtip

# From Pettingill Equation (1), and value for Θ_tip in Table 1.
Θ_tip = 6.9 * pi/180
twist = Θ_tip ./ (r_Rtip)

# Need some aerodynamic quantities.
# Got these using CCBlade.jl: see `AcousticAnalogies.jl/test/gen_bpmjl_data/itr_with_bpmjl.jl`.
data_ccblade = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure22b.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_ccblade["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_ccblade["U"]

# In the text describing Figure 22, "For these predictions, the trip flag was set to “tripped”, due to the rough surface quality of the blade."
bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

# But we're going to use the untripped boundary layer for the LBL-VS noise.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# In the Figure 22 caption, "for these predictions, bluntness thickness H was set to 0.8 mm and trailing edge angle Ψ was set to 16 degrees."
h = 0.8e-3  # meters
Psi = 16*pi/180  # radians

# We'll run for 1 blade pass, 20 time steps per blade pass.
num_blade_pass = 1
num_src_times_blade_pass = 20
bpp = 1/(B/(2*pi)*omega)  # 1/(B blade_passes/rev * 1 rev / (2*pi rad) * omega rad/s)
period_src = num_blade_pass*bpp
num_src_times = num_src_times_blade_pass * num_blade_pass
t0 = 0.0
dt = period_src/num_src_times
src_times = t0 .+ (0:num_src_times-1).*dt

# I don't see any discussion for what type of tip was used for the tip vortex noise.
# The flat tip seems to match the PAS+ROTONET+BARC predictions well.
blade_tip = AcousticAnalogies.FlatTip()

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

# Reshape the inputs to the source element constructors so that everything will line up with (num_times, num_radial, num_blades).
θs_rs = reshape(θs, 1, 1, :)
radii_rs = reshape(radii, 1, :, 1)
dradii_rs = reshape(dradii, 1, :, 1)
# chord_rs = reshape(chord, 1, :, 1)
twist_rs = reshape(twist, 1, :, 1)
# hs_rs = reshape(hs, 1, :, 1)
# Psis_rs = reshape(Psis, 1, :, 1)
Us_rs = reshape(U, 1, :, 1)
alphas_rs = reshape(alpha, 1, :, 1)
# bls_rs = reshape(bls, 1, :, 1)

# Separate things into tip and no-tip.
radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
# chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
# hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
# Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]

radii_rs_with_tip = @view radii_rs[:, end:end, :]
dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
# chord_rs_with_tip = @view chord_rs[:, end:end, :]
twist_rs_with_tip = @view twist_rs[:, end:end, :]
# hs_rs_with_tip = @view hs_rs[:, end:end, :]
# Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
Us_rs_with_tip = @view Us_rs[:, end:end, :]
alphas_rs_with_tip = @view alphas_rs[:, end:end, :]

positive_x_rotation = true
ses_no_tip = CombinedNoTipBroadbandSourceElement.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, Ref(bl), positive_x_rotation) .|> trans
ses_with_tip = CombinedWithTipBroadbandSourceElement.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, Ref(bl), Ref(blade_tip), positive_x_rotation) .|> trans

# It's more convinient to cat all the sources together.
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Do the LBLVS prediction with the untripped boundary layer.
ses_lblvs = LBLVSSourceElement.(asound, nu, radii_rs, θs_rs, dradii_rs, chord, twist_rs, Us_rs, alphas_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# The predictions in Figure 22b appear to be on 1/3 octave, ranging from about 200 Hz to 60,000 Hz.
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
# pbs_lblvss = AcousticAnalogies.pbs_lblvs.(bpm_outs)
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
data_pressure_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-TBL-TE-pressure-2.csv"), ',')
freq_pressure_barc = data_pressure_barc[:, 1]
spl_pressure_barc = data_pressure_barc[:, 2]

data_suction_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-TBL-TE-suction-2.csv"), ',')
freq_suction_barc = data_suction_barc[:, 1]
spl_suction_barc = data_suction_barc[:, 2]

data_separation_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-separation-2.csv"), ',')
freq_separation_barc = data_separation_barc[:, 1]
spl_separation_barc = data_separation_barc[:, 2]

# data_lblvs_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-LBLVS.csv"), ',')
# freq_lblvs_barc = data_lblvs_barc[:, 1]
# spl_lblvs_barc = data_lblvs_barc[:, 2]

data_teb_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-BVS.csv"), ',')
freq_teb_barc = data_teb_barc[:, 1]
spl_teb_barc = data_teb_barc[:, 2]

data_tip_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure22b-tip_vortex_shedding.csv"), ',')
freq_tip_barc = data_tip_barc[:, 1]
spl_tip_barc = data_tip_barc[:, 2]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_pressure_barc, spl_pressure_barc, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_suction_barc, spl_suction_barc, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_separation_barc, spl_separation_barc, color=:yellow, marker=:diamond)
# s_lblvs = scatter!(ax1, freq_lblvs_barc, spl_lblvs_barc, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_teb_barc, spl_teb_barc, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_tip_barc, spl_tip_barc, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
# l_lblvs = lines!(ax1, freqs_obs, spl_lblvs, color=:purple)
l_teb = lines!(ax1, freqs_obs, spl_teb, color=:green)
l_tip = lines!(ax1, freqs_obs, spl_tip, color=:cyan)

xlims!(ax1, 2e2, 6e4)
ylims!(ax1, 10.0, 70.0)

leg = Legend(fig[1, 1], [
        [s_pressure, l_pressure],
        [s_suction, l_suction],
        [s_separation, l_alpha],
        # [s_lblvs, l_lblvs],
        [s_blunt, l_teb],
        [s_tip, l_tip],
    ],
    [
         "TBLTE-Pressure",
         "TBLTE-Suction",
         "Separation",
         # "LBLVS",
         "BVS",
         "Tip",
    ]; orientation=:horizontal, tellwidth=false, tellheight=true, nbanks=2)

text!(ax1, 210, 62; text="markers: CCBlade.jl+BPM.jl\nlines: CCBlade.jl+AcousticAnalogies.jl")

save("figure22b-spl-barc.png", fig)
```
![](figure22b-spl-barc.png)

### Figure 23c
```@example figure23c
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: readdlm
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, ConstantVelocityTransformation
using FileIO: load
using GLMakie
using StaticArrays: @SVector

# Pettingill et al., "Acoustic And Performance Characteristics of an Ideally Twisted Rotor in Hover", 2021
# Parameters from Table 1
B = 4  # number of blades
Rtip = 0.1588 # meters
chord = 0.2*Rtip

# Standard day:
Tamb = 15 + 273.15 # 15°C in Kelvin
pamb = 101325.0  # Pa
R = 287.052874 # J/(kg*K)
rho = pamb/(R*Tamb)
asound = sqrt(1.4*R*Tamb)
# Dynamic and kinematic viscosity
mu = rho*1.4502e-5
nu = mu/rho

# This is a hover case, so the freestream velocity should be zero.
# CCBlade.jl will run with a zero freestream, but I've found that it compares a bit better with experiment if I give it a small non-zero value.
Vinf = 0.001*asound

# Figure 23 caption says Ω_c = 5510 RPM.
rpm = 5510.0
omega = rpm * (2*pi/60)

# Get "cell-centered" radial locations, and also the radial spacing.
num_radial = 50
r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
radii = r_Rtip .* Rtip
dradii = (r_Rtip_[2:end] .- r_Rtip_[1:end-1]) .* Rtip
Rhub = r_Rtip_[1]*Rtip

# From Pettingill Equation (1), and value for Θ_tip in Table 1.
Θ_tip = 6.9 * pi/180
twist = Θ_tip ./ (r_Rtip)

# Need some aerodynamic quantities.
# Got these using CCBlade.jl: see `AcousticAnalogies.jl/test/gen_bpmjl_data/itr_with_bpmjl.jl`.
data_ccblade = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure23c.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_ccblade["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_ccblade["U"]

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

# We'll run for 1 blade pass, 20 time steps per blade pass.
num_blade_pass = 1
num_src_times_blade_pass = 20
bpp = 1/(B/(2*pi)*omega)  # 1/(B blade_passes/rev * 1 rev / (2*pi rad) * omega rad/s)
period_src = num_blade_pass*bpp
num_src_times = num_src_times_blade_pass * num_blade_pass
t0 = 0.0
dt = period_src/num_src_times
src_times = t0 .+ (0:num_src_times-1).*dt

# I don't see any discussion for what type of tip was used for the tip vortex noise.
# The flat tip seems to match the PAS+ROTONET+BARC predictions well.
blade_tip = AcousticAnalogies.FlatTip()

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

# Reshape the inputs to the source element constructors so that everything will line up with (num_times, num_radial, num_blades).
θs_rs = reshape(θs, 1, 1, :)
radii_rs = reshape(radii, 1, :, 1)
dradii_rs = reshape(dradii, 1, :, 1)
# chord_rs = reshape(chord, 1, :, 1)
twist_rs = reshape(twist, 1, :, 1)
# hs_rs = reshape(hs, 1, :, 1)
# Psis_rs = reshape(Psis, 1, :, 1)
Us_rs = reshape(U, 1, :, 1)
alphas_rs = reshape(alpha, 1, :, 1)
bls_rs = reshape(bls, 1, :, 1)

# Separate things into tip and no-tip.
radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
# chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
# hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
# Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]
bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

radii_rs_with_tip = @view radii_rs[:, end:end, :]
dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
# chord_rs_with_tip = @view chord_rs[:, end:end, :]
twist_rs_with_tip = @view twist_rs[:, end:end, :]
# hs_rs_with_tip = @view hs_rs[:, end:end, :]
# Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
Us_rs_with_tip = @view Us_rs[:, end:end, :]
alphas_rs_with_tip = @view alphas_rs[:, end:end, :]
bls_rs_with_tip = @view bls_rs[:, end:end, :]

positive_x_rotation = true
ses_no_tip = CombinedNoTipBroadbandSourceElement.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans
ses_with_tip = CombinedWithTipBroadbandSourceElement.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

# It's more convinient to cat all the sources together.
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS stuff separately.
# Grab the parts of the inputs that correspond to the low Reynolds number stations.
radii_lblvs = @view radii[mask_low_Re_c]
dradii_lblvs = @view dradii[mask_low_Re_c]
# chord_lblvs = @view chord[mask_low_Re_c]
twist_lblvs = @view twist[mask_low_Re_c]
# hs_lblvs = @view hs[mask_low_Re_c]
# Psis_lblvs = @view Psis[mask_low_Re_c]
Us_lblvs = @view U[mask_low_Re_c]
alphas_lblvs = @view alpha[mask_low_Re_c]

# And do the reshaping.
radii_lblvs_rs = reshape(radii_lblvs, 1, :, 1)
dradii_lblvs_rs = reshape(dradii_lblvs, 1, :, 1)
# chord_lblvs_rs = reshape(chord_lblvs, 1, :, 1)
twist_lblvs_rs = reshape(twist_lblvs, 1, :, 1)
# hs_lblvs_rs = reshape(hs_lblvs, 1, :, 1)
# Psis_lblvs_rs = reshape(Psis_lblvs, 1, :, 1)
Us_lblvs_rs = reshape(Us_lblvs, 1, :, 1)
alphas_lblvs_rs = reshape(alphas_lblvs, 1, :, 1)

# Now we can create the source elements.
ses_lblvs = LBLVSSourceElement.(asound, nu, radii_lblvs_rs, θs_rs, dradii_lblvs_rs, chord, twist_lblvs_rs, Us_lblvs_rs, alphas_lblvs_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# Now we can create the source elements.
ses_lblvs = LBLVSSourceElement.(asound, nu, radii_lblvs_rs, θs_rs, dradii_lblvs_rs, chord, twist_lblvs_rs, Us_lblvs_rs, alphas_lblvs_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

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
data_pressure_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-TBL-TE-pressure.csv"), ',')
freq_pressure_barc = data_pressure_barc[:, 1]
spl_pressure_barc = data_pressure_barc[:, 2]

data_suction_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-TBL-TE-suction.csv"), ',')
freq_suction_barc = data_suction_barc[:, 1]
spl_suction_barc = data_suction_barc[:, 2]

data_separation_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-separation.csv"), ',')
freq_separation_barc = data_separation_barc[:, 1]
spl_separation_barc = data_separation_barc[:, 2]

data_lblvs_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-LBLVS.csv"), ',')
freq_lblvs_barc = data_lblvs_barc[:, 1]
spl_lblvs_barc = data_lblvs_barc[:, 2]

data_teb_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-BVS.csv"), ',')
freq_teb_barc = data_teb_barc[:, 1]
spl_teb_barc = data_teb_barc[:, 2]

data_tip_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure23c-tip_vortex_shedding.csv"), ',')
freq_tip_barc = data_tip_barc[:, 1]
spl_tip_barc = data_tip_barc[:, 2]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_pressure_barc, spl_pressure_barc, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_suction_barc, spl_suction_barc, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_separation_barc, spl_separation_barc, color=:yellow, marker=:diamond)
s_lblvs = scatter!(ax1, freq_lblvs_barc, spl_lblvs_barc, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_teb_barc, spl_teb_barc, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_tip_barc, spl_tip_barc, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
l_lblvs = lines!(ax1, freqs_obs, spl_lblvs, color=:purple)
l_teb = lines!(ax1, freqs_obs, spl_teb, color=:green)
l_tip = lines!(ax1, freqs_obs, spl_tip, color=:cyan)

xlims!(ax1, 2e2, 6e4)
ylims!(ax1, 10.0, 70.0)

leg = Legend(fig[1, 1], [
        [s_pressure, l_pressure],
        [s_suction, l_suction],
        [s_separation, l_alpha],
        [s_lblvs, l_lblvs],
        [s_blunt, l_teb],
        [s_tip, l_tip],
    ],
    [
         "TBLTE-Pressure",
         "TBLTE-Suction",
         "Separation",
         "LBLVS",
         "BVS",
         "Tip",
    ]; orientation=:horizontal, tellwidth=false, tellheight=true, nbanks=2)

text!(ax1, 210, 62; text="markers: CCBlade.jl+BPM.jl\nlines: CCBlade.jl+AcousticAnalogies.jl")

save("figure23c-spl-barc.png", fig)
```
![](figure23c-spl-barc.png)


### Figure 24b
```@example figure24b
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using DelimitedFiles: readdlm
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, ConstantVelocityTransformation
using FileIO: load
using GLMakie
using StaticArrays: @SVector

# Pettingill et al., "Acoustic And Performance Characteristics of an Ideally Twisted Rotor in Hover", 2021
# Parameters from Table 1
B = 4  # number of blades
Rtip = 0.1588 # meters
chord = 0.2*Rtip

# Standard day:
Tamb = 15 + 273.15 # 15°C in Kelvin
pamb = 101325.0  # Pa
R = 287.052874 # J/(kg*K)
rho = pamb/(R*Tamb)
asound = sqrt(1.4*R*Tamb)
# Dynamic and kinematic viscosity
mu = rho*1.4502e-5
nu = mu/rho

# This is a hover case, so the freestream velocity should be zero.
# CCBlade.jl will run with a zero freestream, but I've found that it compares a bit better with experiment if I give it a small non-zero value.
Vinf = 0.001*asound

# Figure 24 caption says Ω_c = 2938 RPM.
rpm = 2938.0
omega = rpm * (2*pi/60)

# Get "cell-centered" radial locations, and also the radial spacing.
num_radial = 50
r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
radii = r_Rtip .* Rtip
dradii = (r_Rtip_[2:end] .- r_Rtip_[1:end-1]) .* Rtip
Rhub = r_Rtip_[1]*Rtip

# From Pettingill Equation (1), and value for Θ_tip in Table 1.
Θ_tip = 6.9 * pi/180
twist = Θ_tip ./ (r_Rtip)

# Need some aerodynamic quantities.
# Got these using CCBlade.jl: see `AcousticAnalogies.jl/test/gen_bpmjl_data/itr_with_bpmjl.jl`.
data_ccblade = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure24b.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_ccblade["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_ccblade["U"]

# In the Figure 24 caption, "for these predictions, bluntness thickness H was set to 0.5 mm and trailing edge angle Ψ was set to 14 degrees."
h = 0.5e-3  # meters
Psi = 14*pi/180  # radians

# We'll run for 1 blade pass, 20 time steps per blade pass.
num_blade_pass = 1
num_src_times_blade_pass = 20

# Get the time levels we'll run.
# First, get the blade passing period.
bpp = 1/(B/(2*pi)*omega)  # 1/(B blade_passes/rev * 1 rev / (2*pi rad) * omega rad/s)

# Now we can get the total period of source time we'll run over.
period_src = num_blade_pass*bpp

# And the number of source times.
num_src_times = num_src_times_blade_pass * num_blade_pass

# We know the total time period and number of source times, so we can get the time step.
dt = period_src/num_src_times

# We'll arbitrarily start at time 0.0 seconds.
t0 = 0.0

# And now we can finally get each source time.
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

# For the boundary layer we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
# First figure out how many of each we'll actually have with the `num_radial = 50` radial stations.
num_untripped = Int(round(0.95*num_radial))
num_tripped = num_radial - num_untripped
# Now create a length-`num_radial` vector of untripped and then tripped boundary layer objects.
bls = vcat(fill(AcousticAnalogies.UntrippedN0012BoundaryLayer(), num_untripped), fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), num_tripped))

# But we'll always use the untripped boundary layer with LBLVS.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# I don't see any discussion for what type of tip was used for the tip vortex noise.
# The flat tip seems to match the PAS+ROTONET+BARC predictions well.
blade_tip = AcousticAnalogies.FlatTip()

# Paper doesn't specify the microphone used for Figure 24, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 24.
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

# Reshape the inputs to the source element constructors so that everything will line up with (num_times, num_radial, num_blades).
θs_rs = reshape(θs, 1, 1, :)
radii_rs = reshape(radii, 1, :, 1)
dradii_rs = reshape(dradii, 1, :, 1)
# chord_rs = reshape(chord, 1, :, 1)
twist_rs = reshape(twist, 1, :, 1)
# hs_rs = reshape(hs, 1, :, 1)
# Psis_rs = reshape(Psis, 1, :, 1)
Us_rs = reshape(U, 1, :, 1)
alphas_rs = reshape(alpha, 1, :, 1)
bls_rs = reshape(bls, 1, :, 1)

# Separate things into tip and no-tip.
radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]
bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

radii_rs_with_tip = @view radii_rs[:, end:end, :]
dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
twist_rs_with_tip = @view twist_rs[:, end:end, :]
Us_rs_with_tip = @view Us_rs[:, end:end, :]
alphas_rs_with_tip = @view alphas_rs[:, end:end, :]
bls_rs_with_tip = @view bls_rs[:, end:end, :]

positive_x_rotation = true
ses_no_tip = CombinedNoTipBroadbandSourceElement.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans
ses_with_tip = CombinedWithTipBroadbandSourceElement.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

# Put the source elements together:
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS with the untripped boundary layer.
ses_lblvs = AcousticAnalogies.LBLVSSourceElement.(asound, nu, radii_rs, θs_rs, dradii_rs, chord, twist_rs, Us_rs, alphas_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# The predictions in Figure 24b appear to be on 1/3 octave, ranging from about 200 Hz to 60,000 Hz.
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

# Now I need to account for the fact that Figure 24b is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
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
data_pressure_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-TBL-TE-pressure.csv"), ',')
freq_pressure_barc = data_pressure_barc[:, 1]
spl_pressure_barc = data_pressure_barc[:, 2]

data_suction_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-TBL-TE-suction.csv"), ',')
freq_suction_barc = data_suction_barc[:, 1]
spl_suction_barc = data_suction_barc[:, 2]

data_separation_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-separation.csv"), ',')
freq_separation_barc = data_separation_barc[:, 1]
spl_separation_barc = data_separation_barc[:, 2]

data_lblvs_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-LBLVS.csv"), ',')
freq_lblvs_barc = data_lblvs_barc[:, 1]
spl_lblvs_barc = data_lblvs_barc[:, 2]

data_teb_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-BVS.csv"), ',')
freq_teb_barc = data_teb_barc[:, 1]
spl_teb_barc = data_teb_barc[:, 2]

data_tip_barc = readdlm(joinpath(@__DIR__, "..", "..", "test", "bpm_data", "pettingill_acoustic_performance_characteristics_of_ideally_twisted_rotor_in_hover_2021", "figure24b-tip_vortex_shedding.csv"), ',')
freq_tip_barc = data_tip_barc[:, 1]
spl_tip_barc = data_tip_barc[:, 2]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_pressure_barc, spl_pressure_barc, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_suction_barc, spl_suction_barc, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_separation_barc, spl_separation_barc, color=:yellow, marker=:diamond)
s_lblvs = scatter!(ax1, freq_lblvs_barc, spl_lblvs_barc, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_teb_barc, spl_teb_barc, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_tip_barc, spl_tip_barc, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
l_lblvs = lines!(ax1, freqs_obs, spl_lblvs, color=:purple)
l_teb = lines!(ax1, freqs_obs, spl_teb, color=:green)
l_tip = lines!(ax1, freqs_obs, spl_tip, color=:cyan)

xlims!(ax1, 2e2, 6e4)
ylims!(ax1, 0.0, 50.0)


leg = Legend(fig[1, 1], [
        [s_pressure, l_pressure],
        [s_suction, l_suction],
        [s_separation, l_alpha],
        [s_lblvs, l_lblvs],
        [s_blunt, l_teb],
        [s_tip, l_tip],
    ],
    [
         "TBLTE-Pressure",
         "TBLTE-Suction",
         "Separation",
         "LBLVS",
         "BVS",
         "Tip",
    ]; orientation=:horizontal, tellwidth=false, tellheight=true, nbanks=2)

text!(ax1, 210, 44; text="markers: CCBlade.jl+BPM.jl\nlines: CCBlade.jl+AcousticAnalogies.jl")

save("figure24b-spl-barc.png", fig)
```
![](figure24b-spl-barc.png)

