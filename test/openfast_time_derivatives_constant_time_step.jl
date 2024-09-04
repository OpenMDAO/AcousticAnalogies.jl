using AcousticAnalogies
using Polynomials: fit
using Test

v = azimuth = omega = pitch = nothing
dradii = nothing
cs_area = cs_area_mid = nothing
axial_loading = circum_loading = nothing
T = 2.5
t0 = 0.3
fn(t, r, b) = 0.2*sin(2*pi/T*t*r*b + 0.1*pi) + 0.3*cos(4*pi/T*t*r*b + 0.2*pi)
fndot(t, r, b) = 0.2*2*pi/T*r*b*cos(2*pi/T*t*r*b + 0.1*pi) - 0.3*4*pi/T*r*b*sin(4*pi/T*t*r*b + 0.2*pi)
fc(t, r, b) = 0.4*sin(2*pi/T*t*r*b + 0.3*pi) + 0.5*cos(4*pi/T*t*r*b + 0.4*pi)
fcdot(t, r, b) = 0.4*2*pi/T*r*b*cos(2*pi/T*t*r*b + 0.3*pi) - 0.5*4*pi/T*r*b*sin(4*pi/T*t*r*b + 0.4*pi)
# r = reshape(range(1.0, 2.0; length=5), 1, :)
radii = range(1.0, 2.0; length=5)
radii_mid = 0.5 .* (@view(radii[1:end-1]) .+ @view(radii[2:end]))
r = reshape(radii_mid, 1, :)
b = reshape(1:3, 1, 1, :)
errs_fn_l2 = Vector{Float64}()
errs_fc_l2 = Vector{Float64}()
dts = Vector{Float64}()
for N in 120:10:150
    dt = T/N
    push!(dts, dt)
    time = t0 .+ (0:(N-1)) .* dt
    dtime_dtau = similar(time)
    axial_loading_mid = fn.(time, r, b)
    circum_loading_mid = fc.(time, r, b)
    axial_loading_mid_dot = similar(axial_loading_mid)
    circum_loading_mid_dot = similar(circum_loading_mid)
    # data = OpenFASTData{SecondOrderFiniteDiff}(time, v, azimuth, omega, pitch, axial_loading_mid, circum_loading_mid)
    data = OpenFASTData{FLOWLinearInterp,SecondOrderFiniteDiff}(
        time, dtime_dtau, v, azimuth, omega, pitch,
        radii, radii_mid, dradii,
        cs_area, cs_area_mid,
        axial_loading, axial_loading_mid, axial_loading_mid_dot,
        circum_loading, circum_loading_mid, circum_loading_mid_dot)
    AcousticAnalogies.calculate_loading_dot!(data)
    @test all(data.dtime_dtau .≈ dt)
    err_fn_l2 = sqrt(sum((data.axial_loading_mid_dot .- fndot.(time, r, b)).^2) / length(data.axial_loading_mid_dot))
    push!(errs_fn_l2, err_fn_l2)
    err_fc_l2 = sqrt(sum((data.circum_loading_mid_dot .- fcdot.(time, r, b)).^2) / length(data.circum_loading_mid_dot))
    push!(errs_fc_l2, err_fc_l2)
end

# err ≈ dt^p
# err2/err1 ≈ (dt2^p)/(dt1^p) = (dt2/dt1)^p
# log(err2/err1) ≈ log((dt2/dt1)^p) = p*log(dt2/dt1)
# p ≈ log(err2/err1)/log(dt2/dt1) ≈ (log(err2) - log(err1))/(log(dt2) - log(dt1))

# Fit a line through the errors on a log-log plot, then check that the slope is 2 (second-order).
l_fn = fit(log.(dts), log.(errs_fn_l2), 1)
@test isapprox(l_fn.coeffs[2], 2, atol=0.02)

# Fit a line through the errors on a log-log plot, then check that the slope is 2 (second-order).
l_fc = fit(log.(dts), log.(errs_fc_l2), 1)
@test isapprox(l_fc.coeffs[2], 2, atol=0.02)
