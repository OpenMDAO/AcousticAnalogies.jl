using AcousticAnalogies
using Test

T = 2.5
t0 = 0.3
R = 2.3
num_times = 64
num_radial = 30
num_blades = 3

cs(r) = 0.1*sin(2*pi*r/R + 0.1*pi) + 0.2*cos(4*pi*r/R + 0.2*pi)
fn(t, r, b) = 0.2*sin(2*pi/T*t*r/R*b + 0.1*pi) + 0.3*cos(4*pi/T*t*r/R*b + 0.2*pi)
fc(t, r, b) = 0.4*sin(2*pi/T*t*r/R*b + 0.3*pi) + 0.5*cos(4*pi/T*t*r/R*b + 0.4*pi)

dt = T/num_times
time = t0 .+ (0:(num_times-1)) .* dt

radii = range(0.2*R, R; length=num_radial)
radii_mid = 0.5 .* (@view(radii[1:end-1]) .+ @view(radii[2:end]))
dradii = nothing

dtime_dtau = v = azimuth = omega = pitch = nothing
axial_loading_mid_dot = circum_loading_mid_dot = nothing

r = reshape(radii, 1, :)
b = reshape(1:num_blades, 1, 1, :)
cs_area = @. cs(radii)
axial_loading = @. fn(time, r, b)
circum_loading = @. fc(time, r, b)
cs_area_mid = similar(cs_area, num_radial-1)
axial_loading_mid = similar(axial_loading, num_times, num_radial-1, num_blades)
circum_loading_mid = similar(circum_loading, num_times, num_radial-1, num_blades)

data = OpenFASTData{FLOWLinearInterp,SecondOrderFiniteDiff}(
    time, dtime_dtau, v, azimuth, omega, pitch,
    radii, radii_mid, dradii,
    cs_area, cs_area_mid,
    axial_loading, axial_loading_mid, axial_loading_mid_dot,
    circum_loading, circum_loading_mid, circum_loading_mid_dot)

interpolate_to_cell_centers!(data)

r_mid = reshape(radii_mid, 1, :)
cs_area_mid_exact = @. cs(radii_mid)
axial_loading_mid_exact = @. fn(time, r_mid, b)
circum_loading_mid_exact = @. fc(time, r_mid, b)

csm_min, csm_max = extrema(cs_area_mid_exact)
csm_err = maximum(abs.((data.cs_area_mid .- cs_area_mid_exact)./(csm_max - csm_min)))
@test csm_err < 0.00664

alm_min = reshape(minimum(axial_loading_mid_exact; dims=2), num_times, 1, num_blades)
alm_max = reshape(maximum(axial_loading_mid_exact; dims=2), num_times, 1, num_blades)
alm_err = maximum(abs.((data.axial_loading_mid .- axial_loading_mid_exact)./(alm_max .- alm_min)))
@test alm_err < 0.069

clm_min = reshape(minimum(circum_loading_mid_exact; dims=2), num_times, 1, num_blades)
clm_max = reshape(maximum(circum_loading_mid_exact; dims=2), num_times, 1, num_blades)
clm_err = maximum(abs.((data.circum_loading_mid .- circum_loading_mid_exact)./(clm_max .- clm_min)))
@test clm_err < 0.063
