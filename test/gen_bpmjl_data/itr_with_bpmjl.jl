module ITRWithBPMJL

using Accessors: Accessors
using AcousticMetrics: AcousticMetrics
using BPM: BPM
using CCBlade: CCBlade
using DelimitedFiles: writedlm
# using JLD2: JLD2
using FileIO: save

function get_airfoil(; af_fname, cr75, Re_exp)
    (info, Re, Mach, alpha, cl, cd) = CCBlade.parsefile(af_fname, false)

    # Extend the angle of attack with the Viterna method.
    (alpha, cl, cd) = CCBlade.viterna(alpha, cl, cd, cr75)
    af = CCBlade.AlphaAF(alpha, cl, cd, info, Re, Mach)

    # Reynolds number correction. The 0.6 factor seems to match the NACA 0012
    # drag data from airfoiltools.com.
    reynolds = CCBlade.SkinFriction(Re, Re_exp)

    # Mach number correction.
    mach = CCBlade.PrandtlGlauert()

    # Rotational stall delay correction. Need some parameters from the CL curve.
    m, alpha0 = CCBlade.linearliftcoeff(af, 1.0, 1.0)  # dummy values for Re and Mach
    # Create the Du Selig and Eggers correction.
    rotation = CCBlade.DuSeligEggers(1.0, 1.0, 1.0, m, alpha0)

    # The usual hub and tip loss correction.
    tip = CCBlade.PrandtlTipHub()

    return af, mach, reynolds, rotation, tip
end

function from_cell_centers_to_interfaces(cc_vals::Vector)
    N = length(cc_vals)

    b = cc_vals[:]
    b[1] = 0.25*cc_vals[1] + 0.25*cc_vals[2]

    A = zeros(N, N)
    for i in 1:N-1
        A[i, i] = 0.5
        A[i+1, i] = 0.5
    end
    A[N, N] = 0.5

    a = A\b

    @assert all(A*a .≈ b)

    interface_vals = zeros(N+1)
    interface_vals[1] = 1.5*cc_vals[1] - 0.5*cc_vals[2]
    interface_vals[2:end] .= a

    @assert all(0.5.*(interface_vals[1:end-1] .+ interface_vals[2:end]) .≈ cc_vals)

    return interface_vals
end

function do_figure22b()
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

    # CCBlade.jl defines pitch/collective as "the thing added to the twist to get the orientation of the chord line," but the paper appears to use "the angle of the chord line at the tip."
    collective0 = 6.9*pi/180
    # The Figure 23 caption says Θ_tip = 7 deg, but I think that "actually" means 6.9°.
    # And yes, collective is 0 here, which is a bit silly, I guess, but is consistent with the definition of twist down below.
    collective = 6.9 .* (pi/180) .- collective0

    mu = rho*1.4502e-5
    Vinf = 0.001*asound

    # Figure 22 caption says Ω_c = 5465 RPM.
    rpm = 5465.0
    omega = rpm * (2*pi/60)

    # Get "cell-centered" radial locations.
    num_radial = 50
    r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
    r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
    radii = r_Rtip .* Rtip

    Rhub = r_Rtip_[1]*Rtip

    # From Pettingill Equation (1), and value for Θ_tip in Table 1.
    Θ_tip = 6.9 * pi/180
    twist = Θ_tip ./ (r_Rtip)

    # NACA 0012 airfoil stuff.
    af_fname = joinpath(@__DIR__, "airfoils", "xf-n0012-il-500000.dat")
    Re_exp = 0.6
    cr75 = chord / Rtip
    af, mach_correction, reynolds_correction, rotation_correction, tip_correction = get_airfoil(; af_fname, cr75, Re_exp)
    tip_correction = nothing

    precone = 0.0
    turbine = false
    rotor = CCBlade.Rotor(Rhub, Rtip, B; precone=precone, turbine=turbine, mach=mach_correction, re=reynolds_correction, rotation=rotation_correction, tip=tip_correction)

    sections = CCBlade.Section.(radii, chord, twist, Ref(af))

    ops = CCBlade.OperatingPoint.(Vinf, omega.*radii, rho, collective, mu, asound)
    outs = CCBlade.solve.(Ref(rotor), sections, ops)

    # BPM.jl uses the M_c = 0.8*M assumption from the BPM report, so modify the W field of each CCBlade.Outputs struct to match that.
    lens = Accessors.@optic _.W
    outs_bpm_Mc = Accessors.set.(outs, Ref(lens), 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))
    @assert all(getproperty.(outs_bpm_Mc, :W) .≈ 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))

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
    x0_obs = [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]

    # Get the components of the observer position in a form that BPM.jl uses.
    ox = x0_obs[1]
    oy = x0_obs[2]
    oz = x0_obs[3]

    # Radial locations.
    rs = getproperty.(sections, :r)

    # BPM.jl docstring says we need the angle of attack in degrees.
    # But it expects them at the interfaces, not cell centers.
    alphas_deg = getproperty.(outs_bpm_Mc, :alpha) .* (180/pi)
    rs_interface = from_cell_centers_to_interfaces(rs)
    alphas_deg_interface = from_cell_centers_to_interfaces(alphas_deg)

    # Do the same thing for the chord.
    # What other inputs do we need?
    # So the distance from the pitch axis to the leading edge is used, apparently, to locate the trailing edge.
    # Let's make sure that's true.
    # Here's the relevant code
    #
    #     # Calculate the trailing edge position relative to the hub
    #     xs = sin(beta)*d - cos(beta)*(c - c1)
    #     zs = cos(beta)*d + sin(beta)*(c - c1)
    #
    # `beta` is an azimuthal angle, and `d` is the radial distance of the blade element relative to the hub.
    # So let's say that's zero.
    # Then, I guess, `xs = (c - c1)` and `zs = r`
    # So it looks like the twist is being ignored, I think.
    # But, at any rate, if I set c1 = c, that effect is ignored.
    chords = getproperty.(sections, :chord)
    chords_interface = from_cell_centers_to_interfaces(chords)
    c1s = chords_interface

    # In the text describing Figure 22, "For these predictions, the trip flag was set to “tripped”, due to the rough surface quality of the blade."
    tripped_flags = true

    # In the Figure 22 caption, "for these predictions, bluntness thickness H was set to 0.8 mm and trailing edge angle Ψ was set to 16 degrees."
    h = 0.8e-3  # meters
    Psi = 16*pi/180  # radians
    hs = fill(h, length(rs_interface))
    Psis = fill(Psi, length(rs_interface))
    Psis_deg = Psis .* 180/pi

    # Also need kinematic viscosity.
    nu = mu/rho

    # Number of azimuthal stations per blade pass, I think.
    num_betas = 20

    # Save the inputs.
    data_inputs = Dict(
        "rho"=>rho, "asound"=>asound, "mu"=>mu,
        "Vinf"=>Vinf, "omega"=>omega,
        "B"=>B,
        "Rhub"=>rotor.Rhub,
        "Rtip"=>rotor.Rtip,
        "radii"=>getproperty.(sections, :r),
        "chord"=>getproperty.(sections, :chord),
        "twist"=>getproperty.(sections, :theta),
        "alpha"=>getproperty.(outs, :alpha),
        "U"=>getproperty.(outs, :W),
        "hs"=>hs[begin:end-1],
        "Psis"=>Psis[begin:end-1],
        "tripped_flags"=>fill(tripped_flags, length(sections)),
        "num_src_times_blade_pass"=>num_betas)
    save("figure22b-inputs.jld2", data_inputs)

    # Now we can start doing the predictions.
    oaspl_pressure, spl_pressure = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=true,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_suction, spl_suction = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=true,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_separation, spl_separation = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=true,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_lblvs, spl_lblvs = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=true,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_blunt, spl_blunt = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=true,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_tip, spl_tip = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega, B, rs_interface, chords_interface, c1s, hs, alphas_deg_interface, Psis_deg, nu, asound;
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=true,
        laminar=false,
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    # header = "BPM.default_f,spl_pressure,spl_suction,spl_separation,spl_lblvs,spl_blunt,spl_tip\n"
    # data = hcat(BPM.default_f, spl_pressure, spl_suction, spl_separation, spl_lblvs, spl_blunt, spl_tip)

    # fname = joinpath(@__DIR__, "figure22b.csv")
    # open(fname, "w") do f
    #     write(f, header)
    #     writedlm(f, data, ',')
    # end

    # Now I'd like to do the narrowband SPL like the Pettingill et al. paper does, instead of 1/3 octave SPL.
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
    # First, I'll confirm that the frequencies BPM.jl are using are the ApproximateThirdOctaveCenterBands.
    # This will give me those center bands:
    cbands_approx3rdcenter = AcousticMetrics.ApproximateThirdOctaveCenterBands(first(BPM.default_f), last(BPM.default_f))
    @assert maximum(abs.(BPM.default_f .- cbands_approx3rdcenter)) < 1e-10

    # So then I need to get the spacing associated with the proportional bands.
    # So get the lower and upper "edges" of the bands.
    freqs_l = AcousticMetrics.lower_bands(cbands_approx3rdcenter)
    freqs_u = AcousticMetrics.upper_bands(cbands_approx3rdcenter)
    # And then the spacing for each band.
    df_pbs = freqs_u .- freqs_l
    # Also need the experimental narrowband spacing, which is 20 Hz.
    df_nb = 20.0
    # So, if spl = 10*log10(msp/pref^2), and I want to multiply the msp by df_nb/df_pbs, then
    #
    # spl_nb = 10*log10((msp*df_nb/df_pbs)/pref^2) = 10*(log10(msp/pref^2) + log10(df_nb/df_pbs))
    # spl_nb = 10*log10(msp/pref^2) + 10*log10(df_nb/df_pbs)
    # spl_nb = spl + 10*log10(df_nb/df_pbs)
    # 
    # That's easy.
    spl_nb_pressure = @. spl_pressure + 10*log10(df_nb/df_pbs)
    spl_nb_suction = @. spl_suction + 10*log10(df_nb/df_pbs)
    spl_nb_separation = @. spl_separation + 10*log10(df_nb/df_pbs)
    spl_nb_lblvs = @. spl_lblvs + 10*log10(df_nb/df_pbs)
    spl_nb_blunt = @. spl_blunt + 10*log10(df_nb/df_pbs)
    spl_nb_tip = @. spl_tip + 10*log10(df_nb/df_pbs)

    header = "BPM.default_f,spl_nb_pressure,spl_nb_suction,spl_nb_separation,spl_nb_lblvs,spl_nb_blunt,spl_nb_tip\n"
    data = hcat(BPM.default_f, spl_nb_pressure, spl_nb_suction, spl_nb_separation, spl_nb_lblvs, spl_nb_blunt, spl_nb_tip)

    fname = joinpath(@__DIR__, "figure22b-spl_nb.csv")
    open(fname, "w") do f
        write(f, header)
        writedlm(f, data, ',')
    end
    return nothing
end

function do_figure23c()
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

    # CCBlade.jl defines pitch/collective as "the thing added to the twist to get the orientation of the chord line," but the paper appears to use "the angle of the chord line at the tip."
    collective0 = 6.9*pi/180
    # The Figure 23 caption says Θ_tip = 7 deg, but I think that "actually" means 6.9°.
    # And yes, collective is 0 here, which is a bit silly, I guess, but is consistent with the definition of twist down below.
    collective = 6.9 .* (pi/180) .- collective0

    mu = rho*1.4502e-5
    Vinf = 0.001*asound

    # Figure 23 caption says Ω_c = 5510 RPM.
    rpm = 5510.0
    omega = rpm * (2*pi/60)

    # Get "cell-centered" radial locations.
    num_radial = 50
    r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
    r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
    radii = r_Rtip .* Rtip

    Rhub = r_Rtip_[1]*Rtip

    # From Pettingill Equation (1), and value for Θ_tip in Table 1.
    Θ_tip = 6.9 * pi/180
    twist = Θ_tip ./ (r_Rtip)

    # NACA 0012 airfoil stuff.
    af_fname = joinpath(@__DIR__, "airfoils", "xf-n0012-il-500000.dat")
    Re_exp = 0.6
    cr75 = chord / Rtip
    af, mach_correction, reynolds_correction, rotation_correction, tip_correction = get_airfoil(; af_fname, cr75, Re_exp)
    tip_correction = nothing

    precone = 0.0
    turbine = false
    rotor = CCBlade.Rotor(Rhub, Rtip, B; precone=precone, turbine=turbine, mach=mach_correction, re=reynolds_correction, rotation=rotation_correction, tip=tip_correction)

    sections = CCBlade.Section.(radii, chord, twist, Ref(af))

    ops = CCBlade.OperatingPoint.(Vinf, omega.*radii, rho, collective, mu, asound)
    outs = CCBlade.solve.(Ref(rotor), sections, ops)

    # BPM.jl uses the M_c = 0.8*M assumption from the BPM report, so modify the W field of each CCBlade.Outputs struct to match that.
    lens = Accessors.@optic _.W
    outs_bpm_Mc = Accessors.set.(outs, Ref(lens), 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))
    @assert all(getproperty.(outs_bpm_Mc, :W) .≈ 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))

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
    x0_obs = [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]

    # Get the components of the observer position in a form that BPM.jl uses.
    ox = x0_obs[1]
    oy = x0_obs[2]
    oz = x0_obs[3]

    # Radial locations.
    rs = getproperty.(sections, :r)

    # BPM.jl docstring says we need the angle of attack in degrees.
    # But it expects them at the interfaces, not cell centers.
    rs_interface = from_cell_centers_to_interfaces(rs)
    alphas_deg = getproperty.(outs_bpm_Mc, :alpha) .* (180/pi)
    alphas_deg_interface = from_cell_centers_to_interfaces(alphas_deg)

    # Do the same thing for the chord.
    # What other inputs do we need?
    # So the distance from the pitch axis to the leading edge is used, apparently, to locate the trailing edge.
    # Let's make sure that's true.
    # Here's the relevant code
    #
    #     # Calculate the trailing edge position relative to the hub
    #     xs = sin(beta)*d - cos(beta)*(c - c1)
    #     zs = cos(beta)*d + sin(beta)*(c - c1)
    #
    # `beta` is an azimuthal angle, and `d` is the radial distance of the blade element relative to the hub.
    # So let's say that's zero.
    # Then, I guess, `xs = (c - c1)` and `zs = r`
    # So it looks like the twist is being ignored, I think.
    # But, at any rate, if I set c1 = c, that effect is ignored.
    chords = getproperty.(sections, :chord)
    chords_interface = from_cell_centers_to_interfaces(chords)
    c1s = chords_interface

    # So, for the boundary layer, we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
    # So what are those section indices?
    num_untripped = Int(round(0.95*num_radial))
    num_tripped = num_radial - num_untripped
    _tripped_flags = vcat(fill(false, num_untripped), fill(true, num_tripped))
    @assert length(_tripped_flags) == num_radial
    # This is working with num_radial, the number of cell-centers, not interfaces.
    # I think this is the right thing to do though, because it looks like the last value of the `trip` argument is not used by BPM.jl.
    # But I'll add an extra `true` just to be consistent.
    tripped_flags = vcat(_tripped_flags, _tripped_flags[end:end])

    # Now, the other trick: need to only include LBLVS noise for elements where the Reynolds number is < 160000.
    # So, we need the Reynolds number for each section.
    Re_c = rho * getproperty(outs_bpm_Mc, :W) .* getproperty.(sections, :chord) / mu
    # So now we just need to decide which radial stations have the low or high Re_c values.
    low_Re_c = 160000
    _lblvs_flags = Re_c .< low_Re_c
    # Again, this is working with num_radial-length arrays, which are dealing with cell-centered quantities.
    # Again, the last value is ignored, but I'll add a value to be consistent.
    lblvs_flags = vcat(_lblvs_flags, _lblvs_flags[end:end])

    # In the Figure 23 caption, "for these predictions, bluntness thickness H was set to 0.5 mm and trailing edge angle Ψ was set to 14 degrees."
    h = 0.5e-3  # meters
    Psi = 14*pi/180  # radians
    hs = fill(h, length(rs_interface))
    Psis = fill(Psi, length(rs_interface))
    Psis_deg = Psis .* 180/pi

    # Also need kinematic viscosity.
    nu = mu/rho

    # Number of azimuthal stations per blade pass, I think.
    num_betas = 20

    # Save the inputs.
    data_inputs = Dict(
        "rho"=>rho, "asound"=>asound, "mu"=>mu,
        "Vinf"=>Vinf, "omega"=>omega,
        "B"=>B,
        "Rhub"=>rotor.Rhub,
        "Rtip"=>rotor.Rtip,
        "radii"=>getproperty.(sections, :r),
        "chord"=>getproperty.(sections, :chord),
        "twist"=>getproperty.(sections, :theta),
        "alpha"=>getproperty.(outs, :alpha),
        "U"=>getproperty.(outs, :W),
        "hs"=>hs[begin:end-1],
        "Psis"=>Psis[begin:end-1],
        "tripped_flags"=>tripped_flags[begin:end-1],
        "lblvs_flags"=>lblvs_flags[begin:end-1],
        "num_src_times_blade_pass"=>num_betas)
    save("figure23c-inputs.jld2", data_inputs)

    # Now we can start doing the predictions.
    oaspl_pressure, spl_pressure = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=true,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_suction, spl_suction = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=true,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_separation, spl_separation = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=true,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_lblvs, spl_lblvs = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=lblvs_flags,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_blunt, spl_blunt = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=true,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_tip, spl_tip = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega, B, rs_interface, chords_interface, c1s, hs, alphas_deg_interface, Psis_deg, nu, asound;
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=true,
        laminar=false,
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    # header = "BPM.default_f,spl_pressure,spl_suction,spl_separation,spl_lblvs,spl_blunt,spl_tip\n"
    # data = hcat(BPM.default_f, spl_pressure, spl_suction, spl_separation, spl_lblvs, spl_blunt, spl_tip)

    # fname = joinpath(@__DIR__, "figure23c.csv")
    # open(fname, "w") do f
    #     write(f, header)
    #     writedlm(f, data, ',')
    # end

    # Now I'd like to do the narrowband SPL like the Pettingill et al. paper does, instead of 1/3 octave SPL.
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
    # First, I'll confirm that the frequencies BPM.jl are using are the ApproximateThirdOctaveCenterBands.
    # This will give me those center bands:
    cbands_approx3rdcenter = AcousticMetrics.ApproximateThirdOctaveCenterBands(first(BPM.default_f), last(BPM.default_f))
    @assert maximum(abs.(BPM.default_f .- cbands_approx3rdcenter)) < 1e-10

    # So then I need to get the spacing associated with the proportional bands.
    # So get the lower and upper "edges" of the bands.
    freqs_l = AcousticMetrics.lower_bands(cbands_approx3rdcenter)
    freqs_u = AcousticMetrics.upper_bands(cbands_approx3rdcenter)
    # And then the spacing for each band.
    df_pbs = freqs_u .- freqs_l
    # Also need the experimental narrowband spacing, which is 20 Hz.
    df_nb = 20.0
    # So, if spl = 10*log10(msp/pref^2), and I want to multiply the msp by df_nb/df_pbs, then
    #
    # spl_nb = 10*log10((msp*df_nb/df_pbs)/pref^2) = 10*(log10(msp/pref^2) + log10(df_nb/df_pbs))
    # spl_nb = 10*log10(msp/pref^2) + 10*log10(df_nb/df_pbs)
    # spl_nb = spl + 10*log10(df_nb/df_pbs)
    # 
    # That's easy.
    spl_nb_pressure = @. spl_pressure + 10*log10(df_nb/df_pbs)
    spl_nb_suction = @. spl_suction + 10*log10(df_nb/df_pbs)
    spl_nb_separation = @. spl_separation + 10*log10(df_nb/df_pbs)
    spl_nb_lblvs = @. spl_lblvs + 10*log10(df_nb/df_pbs)
    spl_nb_blunt = @. spl_blunt + 10*log10(df_nb/df_pbs)
    spl_nb_tip = @. spl_tip + 10*log10(df_nb/df_pbs)

    header = "BPM.default_f,spl_nb_pressure,spl_nb_suction,spl_nb_separation,spl_nb_lblvs,spl_nb_blunt,spl_nb_tip\n"
    data = hcat(BPM.default_f, spl_nb_pressure, spl_nb_suction, spl_nb_separation, spl_nb_lblvs, spl_nb_blunt, spl_nb_tip)

    fname = joinpath(@__DIR__, "figure23c-spl_nb.csv")
    open(fname, "w") do f
        write(f, header)
        writedlm(f, data, ',')
    end

    return nothing
end

function do_figure24b()
    # Pettingill et al., "Acoustic And Performance Characteristics of an Ideally Twisted Rotor in Hover", 2021
    # Parameters from Table 1
    B = 4  # number of blades
    Rtip = 0.1588 # meters
    chord = 0.2*Rtip
    # From the first paragraph on page 3, the rotor was designed to produce 11.12 N of thrust.
    # Using that and the value of the design thrust coefficient in Table 1 to get the density and speed of sound.
    # omega_target = 5500.0 * (2*pi/60)
    # thrust_target = 11.12 # Newtons
    # CT_target = 0.0137
    # Mtip_target = 0.27
    # asound = omega_target*Rtip/Mtip_target
    # rho = thrust_target/(CT_target * (pi*Rtip^2) * (omega_target*Rtip)^2)
    # Standard day:
    Tamb = 15 + 273.15 # 15°C in Kelvin
    pamb = 101325.0  # Pa
    R = 287.052874 # J/(kg*K)
    rho = pamb/(R*Tamb)
    asound = sqrt(1.4*R*Tamb)

    # CCBlade.jl defines pitch/collective as "the thing added to the twist to get the orientation of the chord line," but the paper appears to use "the angle of the chord line at the tip."
    collective0 = 6.9*pi/180
    # The Figure 24 caption says Θ_tip = 7 deg, but I think that "actually" means 6.9°.
    # And yes, collective is 0 here, which is a bit silly, I guess, but is consistent with the definition of twist down below.
    collective = 6.9 .* (pi/180) .- collective0

    mu = rho*1.4502e-5
    Vinf = 0.001*asound

    # Figure 24 caption says Ω_c = 2938 RPM.
    rpm = 2938.0
    omega = rpm * (2*pi/60)

    # Get "cell-centered" radial locations.
    num_radial = 50
    r_Rtip_ = range(0.2, 1.0; length=num_radial+1)
    r_Rtip = 0.5 .* (r_Rtip_[2:end] .+ r_Rtip_[1:end-1])
    radii = r_Rtip .* Rtip

    Rhub = r_Rtip_[1]*Rtip

    # From Pettingill Equation (1), and value for Θ_tip in Table 1.
    Θ_tip = 6.9 * pi/180
    twist = Θ_tip ./ (r_Rtip)

    # NACA 0012 airfoil stuff.
    af_fname = joinpath(@__DIR__, "airfoils", "xf-n0012-il-500000.dat")
    Re_exp = 0.6
    cr75 = chord / Rtip
    af, mach_correction, reynolds_correction, rotation_correction, tip_correction = get_airfoil(; af_fname, cr75, Re_exp)
    tip_correction = nothing

    precone = 0.0
    turbine = false
    rotor = CCBlade.Rotor(Rhub, Rtip, B; precone=precone, turbine=turbine, mach=mach_correction, re=reynolds_correction, rotation=rotation_correction, tip=tip_correction)

    sections = CCBlade.Section.(radii, chord, twist, Ref(af))

    ops = CCBlade.OperatingPoint.(Vinf, omega.*radii, rho, collective, mu, asound)
    outs = CCBlade.solve.(Ref(rotor), sections, ops)

    # BPM.jl uses the M_c = 0.8*M assumption from the BPM report, so modify the W field of each CCBlade.Outputs struct to match that.
    lens = Accessors.@optic _.W
    outs_bpm_Mc = Accessors.set.(outs, Ref(lens), 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))
    @assert all(getproperty.(outs_bpm_Mc, :W) .≈ 0.8.*sqrt.(getproperty.(ops, :Vx).^2 .+ getproperty.(ops, :Vy).^2))

    # Paper doesn't specify the microphone used for Figure 24, but earlier at the beginning of "C. Noise Characteristics and Trends" there is this:
    #   > For the purposes of this paper, presented acoustic spectra will correspond to an observer located −35° below the plane of the rotor (microphone 5).
    # So I'll just assume that holds for Figure 24.
    # The observer (microphone 5) is 35 deg behind/downstream of the rotor rotation plane.
    r_obs = 2.27 # meters
    theta_obs = -35*pi/180
    # So, the docstring for BPM.jl says that `V` argument is the wind velocity in the y direction.
    # So I guess we should assume that the blades are rotating about the y axis.
    # And if the freestream velocity is in the positive y axis, then, from the perspective of the fluid, the blades are translating in the negative y direction.
    # And I want the observer to be downstream/behind the blades, so that would mean they would have a positive y position.
    # So I want to rotate the observer around the positive x axis, so I'm going to switch the sign on `theta_obs`.
    t0_obs = 0.0
    x0_obs = [0.0, r_obs*sin(-theta_obs), r_obs*cos(-theta_obs)]

    # Get the components of the observer position in a form that BPM.jl uses.
    ox = x0_obs[1]
    oy = x0_obs[2]
    oz = x0_obs[3]

    # Radial locations.
    rs = getproperty.(sections, :r)

    # BPM.jl docstring says we need the angle of attack in degrees.
    # But it expects them at the interfaces, not cell centers.
    alphas_deg = getproperty.(outs_bpm_Mc, :alpha) .* (180/pi)
    rs_interface = from_cell_centers_to_interfaces(rs)
    alphas_deg_interface = from_cell_centers_to_interfaces(alphas_deg)

    # Do the same thing for the chord.
    # What other inputs do we need?
    # So the distance from the pitch axis to the leading edge is used, apparently, to locate the trailing edge.
    # Let's make sure that's true.
    # Here's the relevant code
    #
    #     # Calculate the trailing edge position relative to the hub
    #     xs = sin(beta)*d - cos(beta)*(c - c1)
    #     zs = cos(beta)*d + sin(beta)*(c - c1)
    #
    # `beta` is an azimuthal angle, and `d` is the radial distance of the blade element relative to the hub.
    # So let's say that's zero.
    # Then, I guess, `xs = (c - c1)` and `zs = r`
    # So it looks like the twist is being ignored, I think.
    # But, at any rate, if I set c1 = c, that effect is ignored.
    chords = getproperty.(sections, :chord)
    chords_interface = from_cell_centers_to_interfaces(chords)
    c1s = chords_interface

    # So, for the boundary layer, we want to use untripped for the 95% of the blade from the hub to almost tip, and then tripped for the last 5% of the blade at the tip.
    # So what are those section indices?
    num_untripped = Int(round(0.95*num_radial))
    num_tripped = num_radial - num_untripped
    _tripped_flags = vcat(fill(false, num_untripped), fill(true, num_tripped))
    @assert length(_tripped_flags) == num_radial
    # This is working with num_radial, the number of cell-centers, not interfaces.
    # I think this is the right thing to do though, because it looks like the last value of the `trip` argument is not used by BPM.jl.
    # But I'll add an extra entry at the end, just to be consistent.
    tripped_flags = vcat(_tripped_flags, _tripped_flags[end:end])

    # In the Figure 24 caption, "for these predictions, bluntness thickness H was set to 0.5 mm and trailing edge angle Ψ was set to 14 degrees."
    h = 0.5e-3  # meters
    Psi = 14*pi/180  # radians
    hs = fill(h, length(rs_interface))
    Psis_deg = fill(Psi*180/pi, length(rs_interface))

    # Also need kinematic viscosity.
    nu = mu/rho

    # Number of azimuthal stations per blade pass, I think.
    num_betas = 20

    # Save the inputs.
    data_inputs = Dict(
        "rho"=>rho, "asound"=>asound, "mu"=>mu,
        "Vinf"=>Vinf, "omega"=>omega,
        "B"=>B,
        "Rhub"=>rotor.Rhub,
        "Rtip"=>rotor.Rtip,
        "radii"=>getproperty.(sections, :r),
        "chord"=>getproperty.(sections, :chord),
        "twist"=>getproperty.(sections, :theta),
        "alpha"=>getproperty.(outs, :alpha),
        "U"=>getproperty.(outs, :W),
        "hs"=>hs[begin:end-1],
        "Psis"=>fill(Psi, length(sections)),
        "tripped_flags"=>tripped_flags[begin:end-1],
        "num_src_times_blade_pass"=>num_betas)
    save("figure24b-inputs.jld2", data_inputs)

    # Now we can start doing the predictions.
    oaspl_pressure, spl_pressure = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=true,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_suction, spl_suction = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=true,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_separation, spl_separation = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=true,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_lblvs, spl_lblvs = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=true,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_blunt, spl_blunt = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega,
        B, rs_interface, chords_interface, c1s, hs,
        alphas_deg_interface, Psis_deg, nu,
        asound;
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        blunt=true,
        weighted=false,
        trip=tripped_flags,
        tip=false,
        laminar=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    oaspl_tip, spl_tip = BPM.sound_pressure_levels(
        ox, oy, oz, Vinf, omega, B, rs_interface, chords_interface, c1s, hs, alphas_deg_interface, Psis_deg, nu, asound;
        blunt=false,
        weighted=false,
        trip=tripped_flags,
        tip=true,
        laminar=false,
        turbulent_pressure=false,
        turbulent_suction=false,
        turbulent_separation=false,
        round=false,
        nbeta=num_betas,
        smooth=false)

    # Now I'd like to do the narrowband SPL like the Pettingill et al. paper does, instead of 1/3 octave SPL.
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
    # First, I'll confirm that the frequencies BPM.jl are using are the ApproximateThirdOctaveCenterBands.
    # This will give me those center bands:
    cbands_approx3rdcenter = AcousticMetrics.ApproximateThirdOctaveCenterBands(first(BPM.default_f), last(BPM.default_f))
    @assert maximum(abs.(BPM.default_f .- cbands_approx3rdcenter)) < 1e-10

    # So then I need to get the spacing associated with the proportional bands.
    # So get the lower and upper "edges" of the bands.
    freqs_l = AcousticMetrics.lower_bands(cbands_approx3rdcenter)
    freqs_u = AcousticMetrics.upper_bands(cbands_approx3rdcenter)
    # And then the spacing for each band.
    df_pbs = freqs_u .- freqs_l
    # Also need the experimental narrowband spacing, which is 20 Hz.
    df_nb = 20.0
    # So, if spl = 10*log10(msp/pref^2), and I want to multiply the msp by df_nb/df_pbs, then
    #
    # spl_nb = 10*log10((msp*df_nb/df_pbs)/pref^2) = 10*(log10(msp/pref^2) + log10(df_nb/df_pbs))
    # spl_nb = 10*log10(msp/pref^2) + 10*log10(df_nb/df_pbs)
    # spl_nb = spl + 10*log10(df_nb/df_pbs)
    # 
    # That's easy.
    spl_nb_pressure = @. spl_pressure + 10*log10(df_nb/df_pbs)
    spl_nb_suction = @. spl_suction + 10*log10(df_nb/df_pbs)
    spl_nb_separation = @. spl_separation + 10*log10(df_nb/df_pbs)
    spl_nb_lblvs = @. spl_lblvs + 10*log10(df_nb/df_pbs)
    spl_nb_blunt = @. spl_blunt + 10*log10(df_nb/df_pbs)
    spl_nb_tip = @. spl_tip + 10*log10(df_nb/df_pbs)

    header = "BPM.default_f,spl_nb_pressure,spl_nb_suction,spl_nb_separation,spl_nb_lblvs,spl_nb_blunt,spl_nb_tip\n"
    data = hcat(BPM.default_f, spl_nb_pressure, spl_nb_suction, spl_nb_separation, spl_nb_lblvs, spl_nb_blunt, spl_nb_tip)

    fname = joinpath(@__DIR__, "figure24b-spl_nb.csv")
    open(fname, "w") do f
        write(f, header)
        writedlm(f, data, ',')
    end

    return nothing
end

function doit()
    do_figure22b()
    do_figure23c()
    do_figure24b()
end

end # module
