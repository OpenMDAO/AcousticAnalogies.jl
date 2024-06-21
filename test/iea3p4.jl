module iea3p4

using AcousticAnalogies
using AcousticMetrics
using ColorSchemes: colorschemes
using GLMakie
using DelimitedFiles
using KinematicCoordinateTransformations
using LinearAlgebra: ×
using StaticArrays
using AcousticMetrics
using Statistics
using Interpolations




## Start of user-defined inputs

# Set number of blades, usually 3 for modern wind turbines
num_blades = 3  # number of blades
# Set the hub radius in m, it is specified in the ElastoDyn main input file of OpenFAST
Rhub = 2.
# Set the blade spanwise grid in m and the corresponding chord, also in m. The two arrays 
# are specified in the AeroDyn15 blade input file
BlSpn = [0.0000e+00, 2.1692e+00, 4.3385e+00, 6.5077e+00, 8.6770e+00, 1.0846e+01, 1.3015e+01, 1.5184e+01, 
    1.7354e+01, 1.9523e+01, 2.1692e+01, 2.3861e+01, 2.6031e+01, 2.8200e+01, 3.0369e+01, 3.2538e+01, 
    3.4708e+01, 3.6877e+01, 3.9046e+01, 4.1215e+01, 4.3385e+01, 4.5554e+01, 4.7723e+01, 
    4.9892e+01, 5.2062e+01, 5.4231e+01, 5.6400e+01, 5.8570e+01, 6.0739e+01, 6.2908e+01]
Chord = [2.600e+00, 2.645e+00, 3.020e+00, 3.437e+00, 3.781e+00, 4.036e+00, 4.201e+00, 
    4.284e+00, 4.288e+00, 4.223e+00, 4.098e+00, 3.923e+00, 3.709e+00, 3.468e+00, 3.220e+00, 
    2.986e+00, 2.770e+00, 2.581e+00, 2.412e+00, 2.266e+00, 2.142e+00, 2.042e+00, 1.964e+00, 
    1.909e+00, 1.870e+00, 1.807e+00, 1.666e+00, 1.387e+00, 9.172e-01, 1.999e-01]

# path to the OpenFAST .out file. The file must contain these channels:
# Time (always available) 
# Wind1VelX from InflowWind
# RotSpeed from ElastoDyn
# Nodal outputs Fxl and Fyl from AeroDyn15
file_path = joinpath(@__DIR__, "gen_test_data", "openfast_data", "IEA-3.4-130-RWT.out")

# Observer location set at the IEC-prescribed location (turbine height on the ground), 
# The microphone is stationary in the global coordinate frame
HH = 110. # m

# Compute the blade locations in radial coordinates, m
RSpn = BlSpn .+ Rhub

x0 = @SVector [HH .+ RSpn[end], 0.0, -HH]

# ambient air density and speed of sound.
rho = 1.225  # kg/m^3
c0 = 340.0  # m/s

## End of user-defined inputs


## Start of computations

# Compute the mid-section locations in radial coordinates, m
radii = 0.5.*(RSpn[begin:end-1] .+ RSpn[begin+1:end])
# Compute the length of each section along blade span, m
dradii = RSpn[begin+1:end] .- RSpn[begin:end-1]
# Compute the blade angles
θs = 2*pi/num_blades.*(0:(num_blades-1))
# Create a linear interpolation object to interpolate chord onto the radial mid-points
itp = LinearInterpolation(RSpn, Chord)
# Perform interpolation
chord = itp(radii)
# Cross-sectional area of each element in m**2. This is taking a bit of a shortcut
cs_area_over_chord_squared = 0.064
cs_area = cs_area_over_chord_squared.*chord.^2

# Code to parse the data from the OpenFAST .out file
# Function to parse a line of data, converting strings to floats
function parse_line(line)
    # Split the line by whitespace and filter out any empty strings
    elements = filter(x -> !isempty(x), split(line))
    # Convert elements to Float64
    return map(x -> parse(Float64, x), elements)
end

# Initialize an empty array to store the data
data = []

# Open the file and read the data, skipping the first 8 lines
open(file_path) do file
    # Skip the first 8 lines (header and description)
    for i in 1:8
        readline(file)
    end

    # Read the rest of the lines and parse them
    for line in eachline(file)
        push!(data, parse_line(line))
    end
end

# Convert the data to an array of arrays (matrix)
data = reduce(hcat, data)
time = data[1, :]
avg_wind_speed = mean(data[2, :])
sim_length_s = time[end] - time[1] # s
@show length(time)

# Reopen the file and read the lines
lines = open(file_path) do f
    readlines(f)
end

# Find the index of the line that contains the column headers
header_index = findfirst(x -> startswith(x, "Time"), lines)

# Extract the headers
headers = split(lines[header_index], '\t')

id_b1_Fx = findfirst(x -> x == "AB1N001Fxl", headers)
id_b2_Fx = findfirst(x -> x == "AB2N001Fxl", headers)
id_b3_Fx = findfirst(x -> x == "AB3N001Fxl", headers)
id_b1_Fy = findfirst(x -> x == "AB1N001Fyl", headers)
id_b2_Fy = findfirst(x -> x == "AB2N001Fyl", headers)
id_b3_Fy = findfirst(x -> x == "AB3N001Fyl", headers)
id_rot_speed = findfirst(x -> x == "RotSpeed", headers)
n_elems = length(radii)
Fx_b1_locs = data[id_b1_Fx:id_b1_Fx+n_elems,:]
Fy_b1_locs = data[id_b1_Fy:id_b1_Fy+n_elems,:]
Fx_b2_locs = data[id_b2_Fx:id_b2_Fx+n_elems,:]
Fy_b2_locs = data[id_b2_Fy:id_b2_Fy+n_elems,:]
Fx_b3_locs = data[id_b3_Fx:id_b3_Fx+n_elems,:]
Fy_b3_locs = data[id_b3_Fy:id_b3_Fy+n_elems,:]

# Reinterpolate onto the mid-sections
Fx_b1 = Array{Float64}(undef, 29, 6001)
Fx_b2 = Array{Float64}(undef, 29, 6001)
Fx_b3 = Array{Float64}(undef, 29, 6001)
Fy_b1 = Array{Float64}(undef, 29, 6001)
Fy_b2 = Array{Float64}(undef, 29, 6001)
Fy_b3 = Array{Float64}(undef, 29, 6001)
for j in axes(Fx_b1_locs, 2)
    itp = LinearInterpolation(RSpn, Fx_b1_locs[:, j], extrapolation_bc=Line())  
    Fx_b1[:, j] = itp(radii)
end
for j in axes(Fx_b2_locs, 2)
    itp = LinearInterpolation(RSpn, Fx_b2_locs[:, j], extrapolation_bc=Line())  
    Fx_b2[:, j] = itp(radii)
end
for j in axes(Fx_b3_locs, 2)
    itp = LinearInterpolation(RSpn, Fx_b3_locs[:, j], extrapolation_bc=Line())  
    Fx_b3[:, j] = itp(radii)
end
for j in axes(Fy_b1_locs, 2)
    itp = LinearInterpolation(RSpn, Fy_b1_locs[:, j], extrapolation_bc=Line())  
    Fy_b1[:, j] = itp(radii)
end
for j in axes(Fy_b2_locs, 2)
    itp = LinearInterpolation(RSpn, Fy_b2_locs[:, j], extrapolation_bc=Line())  
    Fy_b2[:, j] = itp(radii)
end
for j in axes(Fy_b3_locs, 2)
    itp = LinearInterpolation(RSpn, Fy_b3_locs[:, j], extrapolation_bc=Line())  
    Fy_b3[:, j] = itp(radii)
end


# Extract the mean rotor speed
omega_rpm = mean(data[id_rot_speed,:])

# Let's plot the unsteady loading 1 of every 500 timesteps
# x-axis is the span position (mid-sections)
# times are indicated by the colorbar on the right of the plot.
@assert size(Fx_b1) == size(Fx_b2) == size(Fx_b3) == size(Fy_b1) == size(Fy_b2) == size(Fy_b3) 
ntimes_loading = size(Fx_b1, 2)
fig = Figure()
ax11 = fig[1, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Fx (N/m)", title="blade 1")
ax21 = fig[2, 1] = Axis(fig, xlabel="Span Position (m)", ylabel="Fy (N/m)")
ax12 = fig[1, 2] = Axis(fig, xlabel="Span Position (m)", ylabel="Fx (N/m)", title="blade 2")
ax22 = fig[2, 2] = Axis(fig, xlabel="Span Position (m)", ylabel="Fy (N/m)")
ax13 = fig[1, 3] = Axis(fig, xlabel="Span Position (m)", ylabel="Fx (N/m)", title="blade 3")
ax23 = fig[2, 3] = Axis(fig, xlabel="Span Position (m)", ylabel="Fy (N/m)")
bpp = 60/omega_rpm/num_blades
colormap = colorschemes[:viridis]
for tidx in 1:500:ntimes_loading
    cidx = (time[tidx] - time[1])/sim_length_s
    l1 = lines!(ax11, radii, Fx_b1[:,tidx], label ="b1", color=colormap[cidx])
    l1 = lines!(ax12, radii, Fx_b2[:,tidx], label ="b2", color=colormap[cidx])
    l1 = lines!(ax13, radii, Fx_b3[:,tidx], label ="b3", color=colormap[cidx])
    l2 = lines!(ax21, radii, Fy_b1[:,tidx], label ="b1", color=colormap[cidx])
    l2 = lines!(ax22, radii, Fy_b2[:,tidx], label ="b2", color=colormap[cidx])
    l2 = lines!(ax23, radii, Fy_b3[:,tidx], label ="b3", color=colormap[cidx])
end

linkxaxes!(ax21, ax11)
linkxaxes!(ax21, ax11)
linkxaxes!(ax12, ax11)
linkxaxes!(ax22, ax11)
linkxaxes!(ax13, ax11)
linkxaxes!(ax23, ax11)
linkyaxes!(ax12, ax11)
linkyaxes!(ax13, ax11)
linkyaxes!(ax22, ax21)
linkyaxes!(ax23, ax21)

hidexdecorations!(ax11, grid=false)
hidexdecorations!(ax12, grid=false)
hidexdecorations!(ax13, grid=false)
hideydecorations!(ax12, grid=false)
hideydecorations!(ax13, grid=false)
hideydecorations!(ax22, grid=false)
hideydecorations!(ax23, grid=false)

save(joinpath(@__DIR__, "gen_test_data", "openfast_data", "Fx_t-all_time.png"), fig)

# To do F1A correctly, we need to put all source elements in a coordinate system that 
# moves with the fluid, i.e. one in which the fluid appears stationary.
# So, to do that, we have the blades translating in the 
# negative x direction at the average horizontal wind speed.
v = -avg_wind_speed  # m/s
omega = omega_rpm * 2*pi/60  # rad/s

# some reshaping, ses[i, j, k] holds the CompactSourceElement at src_time[i], radii[j], and blade number k
θs = reshape(θs, 1, 1, :)
radii = reshape(radii, 1, :, 1)
dradii = reshape(dradii, 1, :, 1)
cs_area = reshape(cs_area, 1, :, 1)
src_times = reshape(time, :, 1, 1)  # This isn't really necessary.
fx = cat(transpose(Fx_b1), transpose(Fx_b2), transpose(Fx_b3), dims=3)
fc = cat(transpose(Fy_b1), transpose(Fy_b2), transpose(Fy_b3), dims=3)

# source elements, with negative Fx
ses = CompactSourceElement.(rho, c0, radii, θs, dradii, cs_area, -fx, 0.0, fc, src_times)

t0 = 0.0  # Time at which the angle between the source and target coordinate systems is equal to offest.
offset = 0.0  # Angular offset between the source and target coordinate systems at t0.
# steady rotation around the x axis
rot_trans = SteadyRotXTransformation(t0, omega, offset)

# orient the rotation axis of the blades as it is the global frame
rot_axis = @SVector [1.0, 0.0, 0.0] # rotation axis aligned along global x-axis 
blade_axis = @SVector [0.0, 0.0, 1.0]  #  blade 1 pointing up, along z-axis 
global_trans = ConstantLinearMap(hcat(rot_axis, blade_axis, rot_axis×blade_axis))

# blade to move with the appropriate forward velocity, and 
# start from the desired location in the global reference frame
y0_hub = @SVector [0.0, 0.0, 0.0]  # Position of the hub at time t0
v0_hub = SVector{3}(v.*rot_axis)   # Constant velocity of the hub in the global reference frame
const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

# combine these three transformations into one, and then use that on the SourceElements
trans = compose.(src_times, Ref(const_vel_trans), compose.(src_times, Ref(global_trans), Ref(rot_trans)))

# trans will perform the three transformations from right to left (rot_trans, global_trans, const_vel_trans)
ses = ses .|> trans

# The ses object now describes how each blade element source is moving through the global reference
#  frame over the time src_time. As it does this, it will emit acoustics that can be sensed by an acoustic observer 
# (a human, or a microphone). The exact "amount" of acoustics the observer will experience depends 
# on the relative location and motion between each source and the observer. 

# This creates an acoustic observer moving with constant velocity v0_hub that is at location `x0` at time `t0`.
obs = ConstVelocityAcousticObserver(t0, x0, v0_hub)

# Now, in order to perform the F1A calculation, 
# we need to know when each acoustic disturbance emitted 
# by the source arrives at the observer. This is referred 
# to an advanced time calculation, and is done this way:
apth = f1a.(ses, Ref(obs))

# We now have a noise prediction for each of the individual source elements in ses at the acoustic observer obs. 
# What we ultimately want is the total noise prediction at obs—we want to add all the acoustic pressures in apth together. 
# But we can't add them directly, yet, since the observer times are not all the same. What we need to do 
# is first interpolate the apth of each source onto a common observer time grid, and then add them up. 
# We'll do this using the AcousticAnalogies.combine function.
period = 2*pi/omega
bpp = period/num_blades  # blade passing period
obs_time_range = sim_length_s/60*omega_rpm*bpp

# Note that we need to be careful to avoid extrapolation in the `combine` calculation.
# That won't happen in this case, since obs_time_range/sim_length_s is 1/3, so the observer time 
# range is much less than the source time range.
# The observer time range is 1/3 of the source time range, and we're using the same number of 
# simulation times, so that means the observer time step is 1/3 that of the source time step.
num_obs_times = length(time)
apth_total = combine(apth, obs_time_range, num_obs_times, 1)
# The loading data is unsteady, so we may need to be careful to window the time history 
# to avoid problems with discontinuities going from the begining/end of the pressure time history..

# We can now have a look at the total acoustic pressure time history at the observer:
fig = Figure()
ax1 = fig[1, 1] = Axis(fig, xlabel="time, s", ylabel="monopole, Pa")
ax2 = fig[2, 1] = Axis(fig, xlabel="time, s", ylabel="dipole, Pa")
ax3 = fig[3, 1] = Axis(fig, xlabel="time, s", ylabel="total, Pa")
l1 = lines!(ax1, time, apth_total.p_m)
l2 = lines!(ax2, time, apth_total.p_d)
l3 = lines!(ax3, time, apth_total.p_m.+apth_total.p_d)
hidexdecorations!(ax1, grid=false)
hidexdecorations!(ax2, grid=false)
save(joinpath(@__DIR__, "gen_test_data", "openfast_data", "openfast-apth_total.png"), fig)

# The plot shows that the monopole/thickness noise is much lower than the dipole/loading noise.
# Wind turbine blades are relatively slender, which would tend to reduce thickness noise.
# Also the observer is downstream of the rotation plane, which is where loading noise is traditionally 
# thought to dominate (monopole/thickness noise is more significant in the rotor rotation plane, usually).

# We now calculate the overall sound pressure level from the acoustic pressure time history.
oaspl_from_apth = AcousticMetrics.OASPL(apth_total)
# Next, we calculate the narrowband spectrum.
nbs = AcousticMetrics.MSPSpectrumAmplitude(apth_total)
# Finally, walculate the overall A-weighted sound pressure level from the narrowband spectrum.
oaspl_from_nbs = AcousticMetrics.OASPL(nbs)
(oaspl_from_apth, oaspl_from_nbs)

# As a last step, we create VTK files that one can visualize in Paraview (or similar software)
name = joinpath(@__DIR__, "gen_test_data", "openfast_data", "vtk", "iea3p4_vtk")
mkpath(dirname(name))
outfiles = AcousticAnalogies.to_paraview_collection(name, ses)


end # module
