```@meta
CurrentModule = AADocs
```
# Software Quality Assurance, Cont.

## Brooks, Pope, and Marcolini Airfoil Self-Noise Tests, Cont.

### Airfoil Self-Noise Predictions

```@example bpm_figure11_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure11-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0
bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure11-a.png", fig)
```
![](19890016302-figure11-a.png)

```@example bpm_figure11_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure11-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
# M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
M = U/340.46
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure11-b.png", fig)
```
![](19890016302-figure11-b.png)

```@example bpm_figure11_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure11-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
# M = 0.116  # Mach number, corresponds to U = 36.6 m/s in BPM report
M = U/340.46
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure11-c.png", fig)
```
![](19890016302-figure11-c.png)

```@example bpm_figure11_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure11-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# At zero angle of attack the pressure and suction side predictions are the same.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
D_h = AcousticAnalogies.Dbar_h(θ_e, Φ_e, M, M_c)
alphastar = 0.0

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 11 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure11-d.png", fig)
```
![](19890016302-figure11-d.png)

```@example bpm_figure12_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-U71.3-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-U71.3-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-U71.3-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 1.5*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 12 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure12-a.png", fig)
```
![](19890016302-figure12-a.png)

```@example bpm_figure12_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-b-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure12-b-TBL-TE-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 36.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 1.5*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 12 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure12-b.png", fig)
```
![](19890016302-figure12-b.png)

```@example bpm_figure26_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure26-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure26-a.png", fig)
```
![](19890016302-figure26-a.png)

```@example bpm_figure26_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure26-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure26-b.png", fig)
```
![](19890016302-figure26-b.png)

```@example bpm_figure26_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure26-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure26-c.png", fig)
```
![](19890016302-figure26-c.png)

```@example bpm_figure26_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure26-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Pressure and suction sides are the same for zero angle of attack.
f_p = f_s
SPL_p = SPL_s

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # Mach number, corresponds to U = 31.7 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 26 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure26-d.png", fig)
```
![](19890016302-figure26-d.png)

```@example bpm_figure28_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-a-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure28-a.png", fig)
```
![](19890016302-figure28-a.png)

```@example bpm_figure28_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-b-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-b-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 55.5  # freestream velocity in m/s
M = 0.163  # Mach number, corresponds to U = 55.5 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-b.png", fig)
```
![](19890016302-figure28-b.png)

```@example bpm_figure28_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-c-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-c-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-c-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 28 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-c.png", fig)
```
![](19890016302-figure28-c.png)

```@example bpm_figure28_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-d-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure28-d-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 6.7*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="figure 28 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure28-d.png", fig)
```
![](19890016302-figure28-d.png)

```@example bpm_figure38_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure38-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure38-d-TBL-TE-pressure.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_p = bpm[:, 1]
# SPL_p = bpm[:, 2]
# 
# fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure38-d-separation.csv")
# bpm = DelimitedFiles.readdlm(fname, ',')
# f_alpha = bpm[:, 1]
# SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 2.54e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="figure 38 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
# lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")
# 
# scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
# lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure38-d.png", fig)
```
![](19890016302-figure38-d.png)

```@example bpm_figure39_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure39-d-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure39-d-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure39-d-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 2.54e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 4.8*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="figure 39 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure39-d.png", fig)
```
![](19890016302-figure39-d.png)

```@example bpm_figure45_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure45-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure45-a-TBL-TE-pressure.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure45-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure45-a-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 30.48e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 1.5*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 45 (a) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker='◇', label="LBL-VS, BPM")
scatterlines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; marker='◇', label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure45-a.png", fig)
```
![](19890016302-figure45-a.png)

```@example bpm_figure48_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure48-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 22.86e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 0.0*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 48 (c) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 20, 60)
axislegend(ax1, position=:rt)
save("19890016302-figure48-c.png", fig)
```
![](19890016302-figure48-c.png)

```@example bpm_figure54_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure54-a-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 15.24e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 2.7*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 54 (a) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 50, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure54-a.png", fig)
```
![](19890016302-figure54-a.png)

```@example bpm_figure59_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure59-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 0.0*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 59 (c) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure59-c.png", fig)
```
![](19890016302-figure59-c.png)

```@example bpm_figure60_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure60-c-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 3.3*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 60 (c) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure60-c.png", fig)
```
![](19890016302-figure60-c.png)

```@example bpm_figure60_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure60-d-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 10.16e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 3.3*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 60 (d) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure60-d.png", fig)
```
![](19890016302-figure60-d.png)

```@example bpm_figure65_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure65-d-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 31.7  # freestream velocity in m/s
M = 0.093  # mach number, corresponds to u = 31.7 m/s in bpm report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 0.0*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 65 (d) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
lines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 50, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure65-d.png", fig)
```
![](19890016302-figure65-d.png)

```@example bpm_figure66_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure66-b-LBL-VS.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_lbl_vs = bpm[:, 1]
SPL_lbl_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 4.2*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_lbl_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_lblvs=true)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 66 (b) - U = $U m/s")

scatter!(ax1, f_lbl_vs, SPL_lbl_vs; marker=:diamond, label="LBL-VS, BPM")
scatterlines!(ax1, f_jl./1e3, SPL_lbl_vs_jl; label="LBL-VS, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure66-b.png", fig)
```
![](19890016302-figure66-b.png)

```@example bpm_figure69_a
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure69-a-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 71.3  # freestream velocity in m/s
M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 15.4*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 69 (a) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 60, 100)
axislegend(ax1, position=:rt)
save("19890016302-figure69-a.png", fig)
```
![](19890016302-figure69-a.png)

```@example bpm_figure69_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# TBL-TE suction and pressure aren't significant sources for this case (deep stall).

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure69-b-separation.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_alpha = bpm[:, 1]
SPL_alpha = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 5.08e-2  # chord in meters
U = 39.6  # freestream velocity in m/s
M = 0.116  # Mach number, corresponds to U = 39.6 m/s in BPM report
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 15.4*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 69 (b) - U = $U m/s")
# scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

# scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_alpha, SPL_alpha; marker='△', label="separation, BPM")
lines!(ax1, f_jl./1e3, SPL_alpha_jl; label="separation, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure69-b.png", fig)
```
![](19890016302-figure69-b.png)

```@example bpm_figure91
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure91-tip.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_tip = bpm[:, 1]
SPL_tip = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 30.48e-2  # span in meters
chord = 15.24e-2  # chord in meters
speedofsound = 340.46
U = 71.3  # freestream velocity in m/s
# M = 0.209  # Mach number, corresponds to U = 71.3 m/s in BPM report
M = U/speedofsound
M_c = 0.8*M
# speedofsound = U/M
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
alphastar = 10.8*pi/180

bl = AcousticAnalogies.UntrippedN0012BoundaryLayer()
blade_tip = AcousticAnalogies.RoundedTip(AcousticAnalogies.BPMTipAlphaCorrection(), 0.0)
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_tip_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tip_vortex=true, blade_tip=blade_tip)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 91")
scatter!(ax1, f_tip, SPL_tip; marker='o', label="Tip, BPM")
lines!(ax1, f_jl./1e3, SPL_tip_jl; label="Tip, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 90)
axislegend(ax1, position=:rt)
save("19890016302-figure91.png", fig)
```
![](19890016302-figure91.png)

```@example bpm_figure98_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-b-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 1.1e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:rt)
save("19890016302-figure98-b.png", fig)
```
![](19890016302-figure98-b.png)

```@example bpm_figure98_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-c-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 1.9e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:lt)
save("19890016302-figure98-c.png", fig)
```
![](19890016302-figure98-c.png)

```@example bpm_figure98_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 98 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-a-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure98-d-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 69.5  # freestream velocity in m/s
M = U/340.46
h = 2.5e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)


fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 98 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 40, 80)
axislegend(ax1, position=:lt)
save("19890016302-figure98-d.png", fig)
```
![](19890016302-figure98-d.png)

```@example bpm_figure99_b
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 1.1e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (b) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-b.png", fig)
```
![](19890016302-figure99-b.png)

```@example bpm_figure99_c
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-c-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 1.9e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (c) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-c.png", fig)
```
![](19890016302-figure99-c.png)

```@example bpm_figure99_d
using AcousticAnalogies: AcousticAnalogies
using AcousticMetrics: ExactThirdOctaveCenterBands
using DelimitedFiles: DelimitedFiles
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

# Figures 99 a-d only differ in trailing edge bluntness, so the other sources are all the same.
# And TBL-TE is the only significant source, other than bluntness.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_s = bpm[:, 1]
SPL_s = bpm[:, 2]

# Suction and pressure are the same for zero angle of attack.
fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-b-TBL-TE-suction.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_p = bpm[:, 1]
SPL_p = bpm[:, 2]

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure99-d-bluntness.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
f_teb_vs = bpm[:, 1]
SPL_teb_vs = bpm[:, 2]

nu = 1.4529e-5  # kinematic viscosity, m^2/s
L = 45.72e-2  # span in meters
chord = 60.96e-2  # chord in meters
U = 38.6  # freestream velocity in m/s
M = U/340.46
h = 2.5e-3  # trailing edge bluntness in meters
Psi = 14*pi/180  # bluntness angle in radians
r_e = 1.22 # radiation distance in meters
θ_e = 90*pi/180 
Φ_e = 90*pi/180
M_c = 0.8*M
alphastar = 0.0*pi/180

bl = AcousticAnalogies.TrippedN0012BoundaryLayer()
f_jl, SPL_s_jl, SPL_p_jl, SPL_alpha_jl, SPL_teb_vs_jl = AcousticAnalogies.calculate_bpm_test(nu, L, chord, U, M, r_e, θ_e, Φ_e, alphastar, bl; do_tebvs=true, h=h, Psi=Psi)

fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="frequency, kHz", ylabel="SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()),
                       title="Figure 99 (d) - U = $U m/s")
scatter!(ax1, f_s, SPL_s; marker='o', label="TBL-TE suction side, BPM")
lines!(ax1, f_jl./1e3, SPL_s_jl; label="TBL-TE suction side, Julia")

scatter!(ax1, f_p, SPL_p; marker='□', label="TBL-TE pressure side, BPM")
lines!(ax1, f_jl./1e3, SPL_p_jl; label="TBL-TE pressure side, Julia")

scatter!(ax1, f_teb_vs, SPL_teb_vs; marker='◺', label="Bluntness, BPM")
lines!(ax1, f_jl./1e3, SPL_teb_vs_jl; label="Bluntness, Julia")

xlims!(ax1, 0.2, 20.0)
ylims!(ax1, 30, 70)
axislegend(ax1, position=:rt)
save("19890016302-figure99-d.png", fig)
```
![](19890016302-figure99-d.png)
