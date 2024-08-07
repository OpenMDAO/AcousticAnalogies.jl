```@meta
CurrentModule = AADocs
```
# Software Quality Assurance, Cont.

## BPM.jl Comparisons for the Pettingill et al. Ideally Twisted Rotor 
See [here](http://dx.doi.org/10.2514/6.2021-1928) or [here](https://ntrs.nasa.gov/citations/20205003328) for details on the Ideally Twisted Rotor.

### Figure 22b
```@example figure22b
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using GLMakie
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: Akima
using StaticArrays: @SVector

# Copied from BPM.jl (would like to add BPM.jl as a dependency if it's registered in General some day).
# tip vortex noise correction data based on "Airfoil Tip Vortex Formation Noise"
const bm_tip_alpha_aspect_data = [2.0,2.67,4.0,6.0,12.0,24.0]
const bm_tip_alpha_aratio_data = [0.54,0.62,0.71,0.79,0.89,0.95]
const bm_tip_alpha_aspect_ratio_correction = Akima(bm_tip_alpha_aspect_data, bm_tip_alpha_aratio_data)

function bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
    # compute tip lift curve slope
    if aspect_ratio < 2.0
        aratio = 0.5*one(aspect_ratio)
    elseif 2.0 <= aspect_ratio <= 24.0
        aratio = bm_tip_alpha_aspect_ratio_correction(aspect_ratio)
    elseif aspect_ratio > 24.0
        aratio = 1.0*one(aspect_ratio)
    end

    return aratio
end

struct BMTipAlphaCorrection{TCorrection} <: AbstractTipAlphaCorrection
    correction::TCorrection

    function BMTipAlphaCorrection(aspect_ratio)
        # correction = BPM._tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        correction = bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        return new{typeof(correction)}(correction)
    end
end

function AcousticAnalogies.tip_vortex_alpha_correction(blade_tip::AbstractBladeTip{<:BMTipAlphaCorrection}, alphatip)
    a0l = AcousticAnalogies.alpha_zerolift(blade_tip)
    correction_factor = AcousticAnalogies.tip_alpha_correction(blade_tip).correction
    return correction_factor * (alphatip - a0l) + a0l
end

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
data_bpmjl = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure22b.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_bpmjl["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_bpmjl["U"]

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

# BPM.jl uses a different tip alpha correction which appears to require the blade aspect ratio.
# Need to find the blade aspect ratio of the ITR to apply the tip vortex angle of attack correction.
# The aspect ratio is defined as the blade tip radius divided by the average chord, but the chord is constant for this case.
aspect_ratio = Rtip / chord

# Now we can create the tip object.
alpha0lift = 0.0
blade_tip = AcousticAnalogies.FlatTip(BMTipAlphaCorrection(aspect_ratio), alpha0lift)

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

# Use the M_c = 0.8*M that BPM report and BPM.jl use.
U = @. 0.8*sqrt(Vinf^2 + (omega*radii)^2)

# In the text describing Figure 22, "For these predictions, the trip flag was set to “tripped”, due to the rough surface quality of the blade."
# So we'll use a tripped boundary layer for all radial stations along the blade.
bl = AcousticAnalogies.TrippedN0012BoundaryLayer()

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing.
# BPM.jl uses the untripped boundary layer properties for the laminar boundary layer-vortex shedding noise source, so do that here too.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# Paper doesn't specify the microphone used for Figure 22, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 22.
# The observer (microphone 5) is 35 deg behind/downstream of the rotor rotation plane.
r_obs = 2.27 # meters
theta_obs = -35*pi/180
# So, the docstring for BPM.jl says that `V` argument is the wind velocity in the y direction.
# So I guess we should assume that the blades are rotating about the y axis.
# And if the freestream velocity is in the positive y axis, then, from the perspective of the fluid, the blades are translating in the negative y direction.
# And I want the observer to be downstream/behind the blades, so that would mean they would have a positive y position.
# So I want to rotate the observer around the positive x axis, so I'm going to switch the sign on `theta_obs`.
t0_obs = 0.0
x0_obs = @SVector [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]
# The observer is moving in the same direction as the blades, which is the negative y axis.
v_obs = @SVector [0.0, -Vinf, 0.0]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0_obs, x0_obs, v_obs)

# Azimuthal offset for each blade.
θs = (0:(B-1)) .* (2*pi/B) .* ifelse(positive_x_rotation, 1, -1)

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
# bls_untripped_rs = reshape(bls_untripped, 1, :, 1)

# Separate things into tip and no-tip.
radii_rs_no_tip = @view radii_rs[:, begin:end-1, :]
dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
# chord_rs_no_tip = @view chord_rs[:, begin:end-1, :]
twist_rs_no_tip = @view twist_rs[:, begin:end-1, :]
# hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
# Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
Us_rs_no_tip = @view Us_rs[:, begin:end-1, :]
alphas_rs_no_tip = @view alphas_rs[:, begin:end-1, :]
# bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]

radii_rs_with_tip = @view radii_rs[:, end:end, :]
dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
# chord_rs_with_tip = @view chord_rs[:, end:end, :]
twist_rs_with_tip = @view twist_rs[:, end:end, :]
# hs_rs_with_tip = @view hs_rs[:, end:end, :]
# Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
Us_rs_with_tip = @view Us_rs[:, end:end, :]
alphas_rs_with_tip = @view alphas_rs[:, end:end, :]
# bls_rs_with_tip = @view bls_rs[:, end:end, :]

direct = AcousticAnalogies.BPMDirectivity
use_UInduction = false
use_Doppler = false
mach_correction = AcousticAnalogies.NoMachCorrection
ses_no_tip = CombinedNoTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, Ref(bl), positive_x_rotation) .|> trans

ses_with_tip = CombinedWithTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, Ref(bl), Ref(blade_tip), positive_x_rotation) .|> trans

# Put the source elements together:
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing.
lblvs_ses = AcousticAnalogies.LBLVSSourceElement{direct,use_UInduction,use_Doppler}.(asound, nu, radii_rs, θs_rs, dradii_rs, chord, twist_rs, Us_rs, alphas_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans
 
# Define the frequencies we'd like to evaluate.
# BPM.jl uses the approximate 1/3rd-octave bands.
freqs_obs = AcousticMetrics.ApproximateThirdOctaveCenterBands(100.0, 40000.0)
freqs_src = freqs_obs

# Now do the noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))
pbs_lblvss = AcousticAnalogies.noise.(lblvs_ses, Ref(obs), Ref(freqs_src))

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
pbs_lblvs = AcousticMetrics.combine(pbs_lblvss, freqs_obs, time_axis)

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
nb_lblvs_untripped = pbs_lblvs .* df_nb ./ df_pbs

# Now I want the SPL, which should just be this:
pref = 20e-6
spl_pressure = 10 .* log10.(nb_pressure./(pref^2))
spl_suction = 10 .* log10.(nb_suction./(pref^2))
spl_alpha = 10 .* log10.(nb_alpha./(pref^2))
spl_teb = 10 .* log10.(nb_teb./(pref^2))
spl_tip = 10 .* log10.(nb_tip./(pref^2))
spl_lblvs_untripped = 10 .* log10.(nb_lblvs_untripped./(pref^2))

# Finally, let's get the BPM.jl predictions for this case, which we've run and saved previously in a JLD2/HDF5 file.
freq_bpmjl = data_bpmjl["freqs"]
spl_pressure_bpmjl = data_bpmjl["spl_nb_pressure"]
spl_suction_bpmjl = data_bpmjl["spl_nb_suction"]
spl_separation_bpmjl = data_bpmjl["spl_nb_separation"]
spl_lblvs_bpmjl = data_bpmjl["spl_nb_lblvs"]
spl_blunt_bpmjl = data_bpmjl["spl_nb_blunt"]
spl_tip_bpmjl = data_bpmjl["spl_nb_tip"]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_bpmjl, spl_pressure_bpmjl, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_bpmjl, spl_suction_bpmjl, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_bpmjl, spl_separation_bpmjl, color=:yellow, marker=:diamond)
s_lblvs = scatter!(ax1, freq_bpmjl, spl_lblvs_bpmjl, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_bpmjl, spl_blunt_bpmjl, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_bpmjl, spl_tip_bpmjl, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
l_lblvs = lines!(ax1, freqs_obs, spl_lblvs_untripped, color=:purple)
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

save("figure22b-spl-bpmjl.png", fig)
```
![](figure22b-spl-bpmjl.png)

### Figure 23c
```@example figure23c
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using GLMakie
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: Akima
using StaticArrays: @SVector

# Copied from BPM.jl (would like to add BPM.jl as a dependency if it's registered in General some day).
# tip vortex noise correction data based on "Airfoil Tip Vortex Formation Noise"
const bm_tip_alpha_aspect_data = [2.0,2.67,4.0,6.0,12.0,24.0]
const bm_tip_alpha_aratio_data = [0.54,0.62,0.71,0.79,0.89,0.95]
const bm_tip_alpha_aspect_ratio_correction = Akima(bm_tip_alpha_aspect_data, bm_tip_alpha_aratio_data)

function bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
    # compute tip lift curve slope
    if aspect_ratio < 2.0
        aratio = 0.5*one(aspect_ratio)
    elseif 2.0 <= aspect_ratio <= 24.0
        aratio = bm_tip_alpha_aspect_ratio_correction(aspect_ratio)
    elseif aspect_ratio > 24.0
        aratio = 1.0*one(aspect_ratio)
    end

    return aratio
end

struct BMTipAlphaCorrection{TCorrection} <: AbstractTipAlphaCorrection
    correction::TCorrection

    function BMTipAlphaCorrection(aspect_ratio)
        # correction = BPM._tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        correction = bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        return new{typeof(correction)}(correction)
    end
end

function AcousticAnalogies.tip_vortex_alpha_correction(blade_tip::AbstractBladeTip{<:BMTipAlphaCorrection}, alphatip)
    a0l = AcousticAnalogies.alpha_zerolift(blade_tip)
    correction_factor = AcousticAnalogies.tip_alpha_correction(blade_tip).correction
    return correction_factor * (alphatip - a0l) + a0l
end

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
data_bpmjl = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure23c.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_bpmjl["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_bpmjl["U"]

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

# BPM.jl uses a different tip alpha correction which appears to require the blade aspect ratio.
# Need to find the blade aspect ratio of the ITR to apply the tip vortex angle of attack correction.
# The aspect ratio is defined as the blade tip radius divided by the average chord, but the chord is constant for this case.
aspect_ratio = Rtip / chord

# Now we can create the tip object.
alpha0lift = 0.0
blade_tip = AcousticAnalogies.FlatTip(BMTipAlphaCorrection(aspect_ratio), alpha0lift)

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

# Use the M_c = 0.8*M that BPM report and BPM.jl use.
U = @. 0.8*sqrt(Vinf^2 + (omega*radii)^2)

# For the boundary layer we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
# First figure out how many of each we'll actually have with the `num_radial = 50` radial stations.
num_untripped = Int(round(0.95*num_radial))
num_tripped = num_radial - num_untripped
# Now create a length-`num_radial` vector of untripped and then tripped boundary layer objects.
bls = vcat(fill(AcousticAnalogies.UntrippedN0012BoundaryLayer(), num_untripped), fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), num_tripped))

# Now, the other trick: for this case, we're only going to include the LBLVS source where the local Reynolds number (with the chord as the length scale) is < 160000.
low_Re_c = 160000.0
# So we need the Reynolds number for each section.
Re_c = U .* chord / nu
# Now we'll get a length-`num_radial` vector of Bools indicating if that criteria is satisfied or not.
lblvs_flags = Re_c .< low_Re_c

# Also we'll always be using the untripped boundary layer for LBLVS, like BPM.jl does.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# Paper doesn't specify the microphone used for Figure 23, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 23.
# The observer (microphone 5) is 35 deg behind/downstream of the rotor rotation plane.
r_obs = 2.27 # meters
theta_obs = -35*pi/180
# So, the docstring for BPM.jl says that `V` argument is the wind velocity in the y direction.
# So I guess we should assume that the blades are rotating about the y axis.
# And if the freestream velocity is in the positive y axis, then, from the perspective of the fluid, the blades are translating in the negative y direction.
# And I want the observer to be downstream/behind the blades, so that would mean they would have a positive y position.
# So I want to rotate the observer around the positive x axis, so I'm going to switch the sign on `theta_obs`.
t0_obs = 0.0
x0_obs = @SVector [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]
# The observer is moving in the same direction as the blades, which is the negative y axis.
v_obs = @SVector [0.0, -Vinf, 0.0]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0_obs, x0_obs, v_obs)

# Azimuthal offset for each blade.
θs = (0:(B-1)) .* (2*pi/B) .* ifelse(positive_x_rotation, 1, -1)

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

direct = AcousticAnalogies.BPMDirectivity
use_UInduction = false
use_Doppler = false
mach_correction = AcousticAnalogies.NoMachCorrection
ses_no_tip = CombinedNoTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans

ses_with_tip = CombinedWithTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

# Put the source elements together:
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing, and only where `lblvs_flags` is true.
# So extract the radial locations where that's true.
radii_lblvs = @view radii[lblvs_flags]
dradii_lblvs = @view dradii[lblvs_flags]
twist_lblvs = @view twist[lblvs_flags]
Us_lblvs = @view U[lblvs_flags]
alphas_lblvs = @view alpha[lblvs_flags]

# Now do the usual reshaping.
radii_lblvs_rs = reshape(radii_lblvs, 1, :, 1)
dradii_lblvs_rs = reshape(dradii_lblvs, 1, :, 1)
twist_lblvs_rs = reshape(twist_lblvs, 1, :, 1)
Us_lblvs_rs = reshape(Us_lblvs, 1, :, 1)
alphas_lblvs_rs = reshape(alphas_lblvs, 1, :, 1)

# Now we can construct the lblvs source elements.
lblvs_ses = AcousticAnalogies.LBLVSSourceElement{direct,use_UInduction,use_Doppler}.(asound, nu, radii_lblvs_rs, θs_rs, dradii_lblvs_rs, chord, twist_lblvs_rs, Us_lblvs_rs, alphas_lblvs_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# Define the frequencies we'd like to evaluate.
# BPM.jl uses the approximate 1/3rd-octave bands.
freqs_obs = AcousticMetrics.ApproximateThirdOctaveCenterBands(100.0, 40000.0)
freqs_src = freqs_obs

# Now do the noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))
pbs_lblvss = AcousticAnalogies.noise.(lblvs_ses, Ref(obs), Ref(freqs_src))

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
pbs_lblvs = AcousticMetrics.combine(pbs_lblvss, freqs_obs, time_axis)

# Now I need to account for the fact that Figure 23c is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
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
nb_lblvs_untripped = pbs_lblvs .* df_nb ./ df_pbs

# Now I want the SPL, which should just be this:
pref = 20e-6
spl_pressure = 10 .* log10.(nb_pressure./(pref^2))
spl_suction = 10 .* log10.(nb_suction./(pref^2))
spl_alpha = 10 .* log10.(nb_alpha./(pref^2))
spl_teb = 10 .* log10.(nb_teb./(pref^2))
spl_tip = 10 .* log10.(nb_tip./(pref^2))
spl_lblvs_untripped = 10 .* log10.(nb_lblvs_untripped./(pref^2))

# Finally, let's get the BPM.jl predictions for this case, which we've run and saved previously in a JLD2/HDF5 file.
freq_bpmjl = data_bpmjl["freqs"]
spl_pressure_bpmjl = data_bpmjl["spl_nb_pressure"]
spl_suction_bpmjl = data_bpmjl["spl_nb_suction"]
spl_separation_bpmjl = data_bpmjl["spl_nb_separation"]
spl_lblvs_bpmjl = data_bpmjl["spl_nb_lblvs"]
spl_blunt_bpmjl = data_bpmjl["spl_nb_blunt"]
spl_tip_bpmjl = data_bpmjl["spl_nb_tip"]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_bpmjl, spl_pressure_bpmjl, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_bpmjl, spl_suction_bpmjl, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_bpmjl, spl_separation_bpmjl, color=:yellow, marker=:diamond)
s_lblvs = scatter!(ax1, freq_bpmjl, spl_lblvs_bpmjl, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_bpmjl, spl_blunt_bpmjl, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_bpmjl, spl_tip_bpmjl, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
l_lblvs = lines!(ax1, freqs_obs, spl_lblvs_untripped, color=:purple)
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

save("figure23c-spl-bpmjl.png", fig)
```
![](figure23c-spl-bpmjl.png)

### Figure 24b
```@example figure24b
using AcousticAnalogies
using AcousticMetrics: AcousticMetrics
using GLMakie
using KinematicCoordinateTransformations: compose, SteadyRotXTransformation, SteadyRotYTransformation, SteadyRotZTransformation, ConstantVelocityTransformation
using FileIO: load
using FLOWMath: Akima
using StaticArrays: @SVector

# Copied from BPM.jl (would like to add BPM.jl as a dependency if it's registered in General some day).
# tip vortex noise correction data based on "Airfoil Tip Vortex Formation Noise"
const bm_tip_alpha_aspect_data = [2.0,2.67,4.0,6.0,12.0,24.0]
const bm_tip_alpha_aratio_data = [0.54,0.62,0.71,0.79,0.89,0.95]
const bm_tip_alpha_aspect_ratio_correction = Akima(bm_tip_alpha_aspect_data, bm_tip_alpha_aratio_data)

function bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
    # compute tip lift curve slope
    if aspect_ratio < 2.0
        aratio = 0.5*one(aspect_ratio)
    elseif 2.0 <= aspect_ratio <= 24.0
        aratio = bm_tip_alpha_aspect_ratio_correction(aspect_ratio)
    elseif aspect_ratio > 24.0
        aratio = 1.0*one(aspect_ratio)
    end

    return aratio
end

struct BMTipAlphaCorrection{TCorrection} <: AbstractTipAlphaCorrection
    correction::TCorrection

    function BMTipAlphaCorrection(aspect_ratio)
        # correction = BPM._tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        correction = bm_tip_vortex_alpha_correction_nonsmooth(aspect_ratio)
        return new{typeof(correction)}(correction)
    end
end

function AcousticAnalogies.tip_vortex_alpha_correction(blade_tip::AbstractBladeTip{<:BMTipAlphaCorrection}, alphatip)
    a0l = AcousticAnalogies.alpha_zerolift(blade_tip)
    correction_factor = AcousticAnalogies.tip_alpha_correction(blade_tip).correction
    return correction_factor * (alphatip - a0l) + a0l
end

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
data_bpmjl = load(joinpath(@__DIR__, "..", "..", "test", "gen_bpmjl_data", "figure24b.jld2"))
# Angle of attack at each radial station, radians.
alpha = data_bpmjl["alpha"]
# Flow speed normal to span at each radial station, m/s.
U = data_bpmjl["U"]

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
# And th number of source times.
num_src_times = num_src_times_blade_pass * num_blade_pass

# We know the total time period and number of source times, so we can get the time step.
dt = period_src/num_src_times

# We'll arbitrarily start at time 0.0 seconds.
t0 = 0.0

# And now we can finally get each source time.
src_times = t0 .+ (0:num_src_times-1).*dt

# BPM.jl uses a different tip alpha correction which appears to require the blade aspect ratio.
# Need to find the blade aspect ratio of the ITR to apply the tip vortex angle of attack correction.
# The aspect ratio is defined as the blade tip radius divided by the average chord, but the chord is constant for this case.
aspect_ratio = Rtip / chord

# Now we can create the tip object.
alpha0lift = 0.0
blade_tip = AcousticAnalogies.FlatTip(BMTipAlphaCorrection(aspect_ratio), alpha0lift)

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

# Use the M_c = 0.8*M that BPM report and BPM.jl use.
U = @. 0.8*sqrt(Vinf^2 + (omega*radii)^2)

# For the boundary layer we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
# First figure out how many of each we'll actually have with the `num_radial = 50` radial stations.
num_untripped = Int(round(0.95*num_radial))
num_tripped = num_radial - num_untripped
# Now create a length-`num_radial` vector of untripped and then tripped boundary layer objects.
bls = vcat(fill(AcousticAnalogies.UntrippedN0012BoundaryLayer(), num_untripped), fill(AcousticAnalogies.TrippedN0012BoundaryLayer(), num_tripped))

# Also we'll always be using the untripped boundary layer for LBLVS, like BPM.jl does.
bl_lblvs = AcousticAnalogies.UntrippedN0012BoundaryLayer()

# Paper doesn't specify the microphone used for Figure 24, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
#   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
# So I'll just assume that holds for Figure 23.
# The observer (microphone 5) is 35 deg behind/downstream of the rotor rotation plane.
r_obs = 2.27 # meters
theta_obs = -35*pi/180
# So, the docstring for BPM.jl says that `V` argument is the wind velocity in the y direction.
# So I guess we should assume that the blades are rotating about the y axis.
# And if the freestream velocity is in the positive y axis, then, from the perspective of the fluid, the blades are translating in the negative y direction.
# And I want the observer to be downstream/behind the blades, so that would mean they would have a positive y position.
# So I want to rotate the observer around the positive x axis, so I'm going to switch the sign on `theta_obs`.
t0_obs = 0.0
x0_obs = @SVector [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]
# The observer is moving in the same direction as the blades, which is the negative y axis.
v_obs = @SVector [0.0, -Vinf, 0.0]
obs = AcousticAnalogies.ConstVelocityAcousticObserver(t0_obs, x0_obs, v_obs)

# Azimuthal offset for each blade.
θs = (0:(B-1)) .* (2*pi/B) .* ifelse(positive_x_rotation, 1, -1)

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

# Use the directivity functions from the BPM report.
direct = AcousticAnalogies.BPMDirectivity

# Don't include induction for the velocity scale U.
use_UInduction = false

# Don't doppler-shift the source frequencies and source time steps to get observer frequencies & timesteps.
use_Doppler = false

# Don't use the Prandtl-Glauret Mach number correction that Brooks & Burley recommend.
mach_correction = AcousticAnalogies.NoMachCorrection

ses_no_tip = CombinedNoTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_no_tip, θs_rs, dradii_rs_no_tip, chord, twist_rs_no_tip, h, Psi, Us_rs_no_tip, alphas_rs_no_tip, src_times, dt, bls_rs_no_tip, positive_x_rotation) .|> trans

ses_with_tip = CombinedWithTipBroadbandSourceElement{direct,use_UInduction,mach_correction,use_Doppler}.(asound, nu, radii_rs_with_tip, θs_rs, dradii_rs_with_tip, chord, twist_rs_with_tip, h, Psi, Us_rs_with_tip, alphas_rs_with_tip, src_times, dt, bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans

# Put the source elements together:
ses = cat(ses_no_tip, ses_with_tip; dims=2)

# Need to do the LBLVS with the untripped boundary layer to match what BPM.jl is doing.
lblvs_ses = AcousticAnalogies.LBLVSSourceElement{direct,use_UInduction,use_Doppler}.(asound, nu, radii_rs, θs_rs, dradii_rs, chord, twist_rs, Us_rs, alphas_rs, src_times, dt, Ref(bl_lblvs), positive_x_rotation) .|> trans

# Define the frequencies we'd like to evaluate.
# BPM.jl uses the approximate 1/3rd-octave bands.
freqs_obs = AcousticMetrics.ApproximateThirdOctaveCenterBands(100.0, 40000.0)
freqs_src = freqs_obs

# Now do the noise prediction.
bpm_outs = AcousticAnalogies.noise.(ses, Ref(obs), Ref(freqs_src))
pbs_lblvss = AcousticAnalogies.noise.(lblvs_ses, Ref(obs), Ref(freqs_src))

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
pbs_lblvs = AcousticMetrics.combine(pbs_lblvss, freqs_obs, time_axis)

# Now I need to account for the fact that Figure 24b is actually comparing to narrowband experimental data with a frequency spacing of 20 Hz.
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
nb_lblvs_untripped = pbs_lblvs .* df_nb ./ df_pbs

# Now I want the SPL, which should just be this:
pref = 20e-6
spl_pressure = 10 .* log10.(nb_pressure./(pref^2))
spl_suction = 10 .* log10.(nb_suction./(pref^2))
spl_alpha = 10 .* log10.(nb_alpha./(pref^2))
spl_teb = 10 .* log10.(nb_teb./(pref^2))
spl_tip = 10 .* log10.(nb_tip./(pref^2))
spl_lblvs_untripped = 10 .* log10.(nb_lblvs_untripped./(pref^2))

# Finally, let's get the BPM.jl predictions for this case, which we've run and saved previously in a JLD2/HDF5 file.
freq_bpmjl = data_bpmjl["freqs"]
spl_pressure_bpmjl = data_bpmjl["spl_nb_pressure"]
spl_suction_bpmjl = data_bpmjl["spl_nb_suction"]
spl_separation_bpmjl = data_bpmjl["spl_nb_separation"]
spl_lblvs_bpmjl = data_bpmjl["spl_nb_lblvs"]
spl_blunt_bpmjl = data_bpmjl["spl_nb_blunt"]
spl_tip_bpmjl = data_bpmjl["spl_nb_tip"]

# Now let's plot.
fig = Figure()
ax1 = fig[2, 1] = Axis(fig, xlabel="frequency, Hz", ylabel="SPL (dB Ref: 20 μPa), Δf = 20 Hz", xscale=log10, xticks=[10^3, 10^4], xminorticksvisible=true, xminorgridvisible=true, xminorticks=IntervalsBetween(9), yticks=10:10:70)#, aspect=3)

s_pressure = scatter!(ax1, freq_bpmjl, spl_pressure_bpmjl, color=:blue, marker=:rtriangle)
s_suction = scatter!(ax1, freq_bpmjl, spl_suction_bpmjl, color=:red, marker=:ltriangle)
s_separation = scatter!(ax1, freq_bpmjl, spl_separation_bpmjl, color=:yellow, marker=:diamond)
s_lblvs = scatter!(ax1, freq_bpmjl, spl_lblvs_bpmjl, color=:purple, marker=:rect)
s_blunt = scatter!(ax1, freq_bpmjl, spl_blunt_bpmjl, color=:green, marker=:star6)
s_tip = scatter!(ax1, freq_bpmjl, spl_tip_bpmjl, color=:cyan, marker=:circle)

l_pressure = lines!(ax1, freqs_obs, spl_pressure, color=:blue)
l_suction = lines!(ax1, freqs_obs, spl_suction, color=:red)
l_alpha = lines!(ax1, freqs_obs, spl_alpha, color=:yellow)
l_lblvs = lines!(ax1, freqs_obs, spl_lblvs_untripped, color=:purple)
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

save("figure24b-spl-bpmjl.png", fig)
```
![](figure24b-spl-bpmjl.png)
