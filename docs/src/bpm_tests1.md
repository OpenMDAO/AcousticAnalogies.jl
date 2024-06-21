```@meta
CurrentModule = AADocs
```
# Software Quality Assurance, Cont.

## Brooks, Pope, and Marcolini Airfoil Self-Noise Tests
The [Brooks, Pope, and Marcolini (BPM) report on airfoil self-noise](https://ntrs.nasa.gov/citations/19890016302) forms the basis of the [Brooks and Burley broadband noise modeling approach](https://doi.org/10.2514/6.2001-2210) that is implemented in AcousticAnalogies.jl.

### Boundary Layer Tests

```@example bpm_bl_thickness
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/10^6", ylabel="δ_0/c",
                       xscale=log10, yscale=log10,
                       xminorticksvisible=true, yminorticksvisible=true,
                       xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()), yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-bl_thickness-tripped.csv")
bpm_tripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_tripped[:, 1]
deltastar0_c = bpm_tripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="tripped, BPM report", color=colors[1])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.bl_thickness_0.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
@show Re_c_1e6_jl deltastar0_c_jl
@show typeof(Re_c_1e6_jl) typeof(deltastar0_c_jl)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="tripped, Julia", color=colors[1])

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-bl_thickness-untripped.csv")
bpm_untripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_untripped[:, 1]
deltastar0_c = bpm_untripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="untripped, BPM report", color=colors[2])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.bl_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="untripped, Julia", color=colors[2])

xlims!(ax1, 0.04, 3)
ylims!(ax1, 0.01, 0.2)
axislegend(ax1)
save("19890016302-figure06-bl_thickness.png", fig)
```
![](19890016302-figure06-bl_thickness.png)

```@example bpm_disp_thickness
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c/10^6", ylabel="δ_0^*/c",
                       xscale=log10, yscale=log10,
                       xminorticksvisible=true, yminorticksvisible=true,
                       xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()), yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-disp_thickness-tripped.csv")
bpm_tripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_tripped[:, 1]
deltastar0_c = bpm_tripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="tripped, BPM report", color=colors[1])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="tripped, Julia", color=colors[1])

# deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)
# lines!(ax1, Re_c_1e6_jl, deltastar0_c_interp; linestyle=:dash, label="tripped, interp", color=colors[1])

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure06-disp_thickness-untripped.csv")
bpm_untripped = DelimitedFiles.readdlm(fname, ',')
Re_c_1e6 = bpm_untripped[:, 1]
deltastar0_c = bpm_untripped[:, 2]
scatter!(ax1, Re_c_1e6, deltastar0_c, markersize=4, label="untripped, BPM report", color=colors[2])

Re_c_1e6_jl = range(minimum(Re_c_1e6), maximum(Re_c_1e6); length=50)
deltastar0_c_jl = AcousticAnalogies.disp_thickness_0.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), Re_c_1e6_jl.*1e6)
lines!(ax1, Re_c_1e6_jl, deltastar0_c_jl, label="untripped, Julia", color=colors[2])

# deltastar0_c_interp = linear(Re_c_1e6, deltastar0_c, Re_c_1e6_jl)
# lines!(ax1, Re_c_1e6_jl, deltastar0_c_interp; linestyle=:dash, label="untripped, interp", color=colors[2])

xlims!(ax1, 0.04, 3)
ylims!(ax1, 0.001, 0.03)
axislegend(ax1)
save("19890016302-figure06-disp_thickness.png", fig)
```
![](19890016302-figure06-disp_thickness.png)

```@example bpm_bl_thickness_tripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ/δ_0",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-bl_thickness-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
delta_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, delta_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
delta_jl = AcousticAnalogies._bl_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, delta_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, delta_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(delta_bpm)
# err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-bl_thickness-suction_side.csv")
bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_suction_side[:, 1]
delta_bpm = bpm_suction_side[:, 2]
scatter!(ax1, alpha_deg, delta_bpm, markersize=4, color=colors[2], label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
delta_jl = AcousticAnalogies._bl_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, delta_jl; color=colors[2], label="suction side, Julia")

# Interpolate:
# delta_bpm_interp = linear(alpha_deg, delta_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, delta_bpm_interp, color=colors[2], linestyle=:dash, label="suction side, interp")
# Check error.
# vmin, vmax = extrema(delta_bpm)
# err = abs.(delta_jl .- delta_bpm_interp)./(vmax - vmin)
# println("suction side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 20)
axislegend(ax1, position=:lt)
save("19890016302-figure07-bl_thickness.png", fig)
```
![](19890016302-figure07-bl_thickness.png)

```@example bpm_disp_thickness_star_tripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ^*/δ_0^*",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._disp_thickness_p.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure07-suction_side.csv")
bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_suction_side[:, 1]
deltastar_bpm = bpm_suction_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, markersize=4, color=colors[2], label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._disp_thickness_s.(Ref(AcousticAnalogies.TrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[2], label="suction side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[2], linestyle=:dash, label="suction side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("suction side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 200)
axislegend(ax1, position=:lt)
save("19890016302-figure07.png", fig)
```
![](19890016302-figure07.png)

```@example bpm_bl_thickness_untripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ/δ_0",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-bl_thickness-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._bl_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-bl_thickness-suction_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[2], markersize=4, label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._bl_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[2], label="suction side, Julia")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 40)
axislegend(ax1, position=:lt)
save("19890016302-figure08-bl_thickness.png", fig)
```
![](19890016302-figure08-bl_thickness.png)

```@example bpm_disp_thickness_star_untripped
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="alpha, deg.", ylabel="δ^*/δ_0^*",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks())
                       )

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-pressure_side.csv")
bpm_pressure_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_pressure_side[:, 1]
deltastar_bpm = bpm_pressure_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, color=colors[1], markersize=4, label="pressure side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._disp_thickness_p.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[1], label="pressure side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[1], linestyle=:dash, label="pressure side, interp")
# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("pressure side error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure08-suction_side.csv")
bpm_suction_side = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm_suction_side[:, 1]
deltastar_bpm = bpm_suction_side[:, 2]
scatter!(ax1, alpha_deg, deltastar_bpm, markersize=4, color=colors[2], label="suction side, BPM report")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
deltastar_jl = AcousticAnalogies._disp_thickness_s.(Ref(AcousticAnalogies.UntrippedN0012BoundaryLayer()), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, deltastar_jl; color=colors[2], label="suction side, Julia")

# Interpolate:
# deltastar_bpm_interp = linear(alpha_deg, deltastar_bpm, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, deltastar_bpm_interp, color=colors[2], linestyle=:dash, label="suction side, interp")

# Check error.
# vmin, vmax = extrema(deltastar_bpm)
# err = abs.(deltastar_jl .- deltastar_bpm_interp)./(vmax - vmin)
# println("suction side error =\n$(err)")

xlims!(ax1, 0, 25)
ylims!(ax1, 0.2, 200)
axislegend(ax1, position=:lt)
save("19890016302-figure08.png", fig)
```
![](19890016302-figure08.png)

### Turbulent Boundary Layer-Trailing Edge Tests

```@example bpm_K_1
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Re_c", ylabel="Peak scaled SPL_1/3, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure77.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
Re_c_bpm = bpm[:, 1]
K_1_bpm = bpm[:, 2]
scatter!(ax1, Re_c_bpm, K_1_bpm, color=colors[1], markersize=8, label="BPM report")

Re_c_jl = range(minimum(Re_c_bpm), maximum(Re_c_bpm); length=50)
K_1_jl = AcousticAnalogies.K_1.(Re_c_jl)
lines!(ax1, Re_c_jl, K_1_jl, color=colors[1], label="Julia")

# Interpolate:
# K_1_interp = linear(Re_c_bpm, K_1_bpm, Re_c_jl)
# lines!(ax1, Re_c_jl, K_1_interp; color=colors[1], linestyle=:dash, label="interp")

# # Check error.
# vmin, vmax = extrema(K_1_bpm)
# err = abs.(K_1_jl .- K_1_interp)./(vmax - vmin)
# println("K_1_error =\n$(err)")

xlims!(ax1, 10^4, 10^7)
ylims!(ax1, 110.0, 150.0)
axislegend(ax1, position=:lt)
save("19890016302-figure77.png", fig)
```
![](19890016302-figure77.png)

```@example bpm_A
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal number ratio, St/St_peak", ylabel="Function A level, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-A_min.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
A = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, A, color=colors[1], markersize=8, label="A_min, BPM report")

St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
A_jl = AcousticAnalogies.A.(St_St_peak_jl, 9.5e4)
lines!(ax1, St_St_peak_jl, A_jl, color=colors[1], label="A_min, Julia")

# # Interpolate:
# A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, A_interp; color=colors[1], linestyle=:dash, label="A_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(A)
# err = abs.(A_jl .- A_interp)./(vmax - vmin)
# println("A_min error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-A_max.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
A = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, A, color=colors[2], markersize=8, label="A_max, BPM report")

St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
A_jl = AcousticAnalogies.A.(St_St_peak_jl, 8.58e5)
lines!(ax1, St_St_peak_jl, A_jl, color=colors[2], label="A_max, Julia")

# # Interpolate:
# A_interp = linear(St_St_peak_bpm, A, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, A_interp; color=colors[2], linestyle=:dash, label="A_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(A)
# err = abs.(A_jl .- A_interp)./(vmax - vmin)
# println("A_max error =\n$(err)")

xlims!(ax1, 0.1, 20)
ylims!(ax1, -20.0, 0.0)
axislegend(ax1, position=:lt)
save("19890016302-figure78-A.png", fig)
```
![](19890016302-figure78-A.png)

```@example bpm_B
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Strouhal number ratio, St/St_peak", ylabel="Function B level, dB",
                       xscale=log10,
                       xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9),
                       xticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-B_min.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
B = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, B, color=colors[1], markersize=8, label="B_min, BPM report")

# St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
St_St_peak_jl = range(0.5, 2; length=50)
B_jl = AcousticAnalogies.B.(St_St_peak_jl, 9.5e4)
lines!(ax1, St_St_peak_jl, B_jl, color=colors[1], label="B_min, Julia")

# # Interpolate:
# B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, B_interp; color=colors[1], linestyle=:dash, label="B_min, interp")
# 
# # Check error.
# vmin, vmax = extrema(B)
# err = abs.(B_jl .- B_interp)./(vmax - vmin)
# println("B_min error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure78-B_max.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
St_St_peak_bpm = bpm[:, 1]
B = bpm[:, 2]
scatter!(ax1, St_St_peak_bpm, B, color=colors[2], markersize=8, label="B_max, BPM report")

# St_St_peak_jl = range(minimum(St_St_peak_bpm), maximum(St_St_peak_bpm); length=50)
St_St_peak_jl = range(0.2, 4; length=50)
B_jl = AcousticAnalogies.B.(St_St_peak_jl, 8.58e5)
lines!(ax1, St_St_peak_jl, B_jl, color=colors[2], label="B_max, Julia")

# # Interpolate:
# B_interp = linear(St_St_peak_bpm, B, St_St_peak_jl)
# lines!(ax1, St_St_peak_jl, B_interp; color=colors[2], linestyle=:dash, label="B_max, interp")
# 
# # Check error.
# vmin, vmax = extrema(B)
# err = abs.(B_jl .- B_interp)./(vmax - vmin)
# println("B_max error =\n$(err)")

xlims!(ax1, 0.1, 20)
ylims!(ax1, -20.0, 0.0)
axislegend(ax1, position=:lt)
save("19890016302-figure78-B.png", fig)
```
![](19890016302-figure78-B.png)

```@example bpm_St_2
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Angle of attack α^*, deg", ylabel="Peak Strouhal number, St_peak",
                       yscale=log10,
                       yminorticksvisible=true,
                       yminorticks=IntervalsBetween(9),
                       yticks=LogTicks(IntegerTicks()))

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure80-M0.093.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
St_2 = bpm[:, 2]
scatter!(ax1, alpha_deg, St_2, color=colors[1], markersize=8, label="St_2 for M = 0.093, BPM")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.093), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, St_2_jl, color=colors[1], label="St_2 for M = 0.093, Julia")

# # Interpolate:
# St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, St_2_interp; color=colors[1], linestyle=:dash, label="St_2 for M = 0.093, interp")
# 
# # Check error.
# vmin, vmax = extrema(St_2)
# err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
# println("St_2 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure80-M0.209.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
St_2 = bpm[:, 2]
scatter!(ax1, alpha_deg, St_2, color=colors[2], markersize=8, label="St_2 for M = 0.209, BPM")

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=50)
St_2_jl = AcousticAnalogies.St_2.(AcousticAnalogies.St_1(0.209), alpha_deg_jl.*pi/180)
lines!(ax1, alpha_deg_jl, St_2_jl, color=colors[2], label="St_2 for M = 0.209, Julia")

# # Interpolate:
# St_2_interp = linear(alpha_deg, St_2, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, St_2_interp; color=colors[2], linestyle=:dash, label="St_2 for M = 0.209, interp")
# 
# # Check error.
# vmin, vmax = extrema(St_2)
# err = abs.(St_2_jl .- St_2_interp)./(vmax - vmin)
# println("St_2 error =\n$(err)")

xlims!(ax1, 0.0, 25.0)
ylims!(ax1, 0.01, 1)
axislegend(ax1, position=:lt)
save("19890016302-figure80.png", fig)
```
![](19890016302-figure80.png)

```@example bpm_K_2_K_1
using AcousticAnalogies: AcousticAnalogies
using ColorSchemes: colorschemes
using DelimitedFiles: DelimitedFiles
# using FLOWMath: linear
using GLMakie

# https://docs.makie.org/stable/examples/blocks/axis/index.html#logticks
struct IntegerTicks end
Makie.get_tickvalues(::IntegerTicks, vmin, vmax) = ceil(Int, vmin) : floor(Int, vmax)

colors = colorschemes[:tab10]
fig = Figure()
ax1 = fig[1, 1] = Axis(fig; xlabel="Angle of attack α_*, deg", ylabel="Extracted scaled levels minus K_1, dB")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.093.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[1], markersize=8, label="M = 0.093, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.093, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[1], label="M = 0.093, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[1], linestyle=:dash, label="M = 0.093, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.116.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[2], markersize=8, label="M = 0.116, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.116, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[2], label="M = 0.116, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[2], linestyle=:dash, label="M = 0.116, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.163.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[3], markersize=8, label="M = 0.163, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.163, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[3], label="M = 0.163, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[3], linestyle=:dash, label="M = 0.163, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

fname = joinpath(@__DIR__, "..", "..", "test", "bpm_data", "brooks_airfoil_self_noise_and_prediction_1989", "19890016302-figure82-M0.209.csv")
bpm = DelimitedFiles.readdlm(fname, ',')
alpha_deg = bpm[:, 1]
K_2_K_1 = bpm[:, 2]
scatter!(ax1, alpha_deg, K_2_K_1, color=colors[4], markersize=8, label="M = 0.209, BPM", marker='o')

alpha_deg_jl = range(minimum(alpha_deg), maximum(alpha_deg); length=200)
K_2_K_1_jl = AcousticAnalogies.K_2.(1e6, 0.209, alpha_deg_jl.*pi/180) .- AcousticAnalogies.K_1(1e6)
lines!(ax1, alpha_deg_jl, K_2_K_1_jl, color=colors[4], label="M = 0.209, Julia")

# # Interpolate:
# K_2_K_1_interp = linear(alpha_deg, K_2_K_1, alpha_deg_jl)
# lines!(ax1, alpha_deg_jl, K_2_K_1_interp; color=colors[4], linestyle=:dash, label="M = 0.209, interp")
# 
# # Check error.
# vmin, vmax = extrema(K_2_K_1)
# err = abs.(K_2_K_1_jl .- K_2_K_1_interp)./(vmax - vmin)
# println("K_2_K_1 error =\n$(err)")

xlims!(ax1, 0.0, 25.0)
ylims!(ax1, -20, 20)
axislegend(ax1, position=:lt)
save("19890016302-figure82.png", fig)
```
![](19890016302-figure82.png)

