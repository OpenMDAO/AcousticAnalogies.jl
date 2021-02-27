module GenANOPP2Data

import ANOPP2, DelimitedFiles
using Printf: @sprintf

include("gen_ccblade_data/constants.jl")

"""
    get_dradii(radii, Rhub, Rtip)

Compute the spacing between blade elements given the radial locations of the
element midpoints in `radii` and the hub and tip radius in `Rhub` and `Rtip`,
respectively.
"""
function get_dradii(radii, Rhub, Rtip)
    # How do I get the radial spacing? Well, for the inner elements, I'll just
    # assume that the interfaces are midway in between the centers.
    r_interface = 0.5.*(radii[1:end-1] .+ radii[2:end])
    # Then just say that the blade begins at Rhub, and ends at Rtip.
    r_interface = vcat([Rhub], r_interface, [Rtip])
    # And now the distance between interfaces is the spacing.
    dradii = r_interface[2:end] .- r_interface[1:end-1]
    return dradii
end

function omega_sweep_with_acoustics(; stationary_observer, theta)
    dradii = get_dradii(radii, Rhub, Rtip)
    cs_area =  area_over_chord_squared .* chord.^2

    rpm = 200.0:200.0:2200.0  # rev/min
    for i in eachindex(rpm)
        omega = rpm[i]*(2*pi/60.0)

        # Get the normal and circumferential loading from the CCBlade output.
        data = DelimitedFiles.readdlm("gen_ccblade_data/ccblade_omega$(@sprintf "%02d" i).csv", ',')
        fn = data[:, 1]
        fc = data[:, 2]

        # Calculate the noise with ANOPP2.
        t_monopole, p_monopole, t_dipole, p_dipole = anopp2_noise(num_blades, v, omega, radii, dradii, chord, fn, fc, stationary_observer, theta)

        # Check that the two times are the same.
        max_diff = maximum(abs.(t_monopole .- t_dipole))
        if max_diff > 1e-10
            msg = """
                time vectors associated with the monopole and dipole sources appear to differ"
                t_monopole = $(t_monopole)
                t_dipole   = $(t_dipole)
            """
            @warn msg
        end

        # Write the data.
        data = hcat(t_monopole, p_monopole, p_dipole)
        if theta ≈ 0.0
            if stationary_observer
                DelimitedFiles.writedlm("anopp2_omega$(@sprintf "%02d" i).csv", data, ',')
            else
                DelimitedFiles.writedlm("anopp2_const_vel_omega$(@sprintf "%02d" i).csv", data, ',')
            end
        else
            theta_int = Int(round(theta*180.0/pi))
            if stationary_observer
                DelimitedFiles.writedlm("anopp2_omega$(@sprintf "%02d" i)_theta$(theta_int).csv", data, ',')
            else
                DelimitedFiles.writedlm("anopp2_const_vel_omega$(@sprintf "%02d" i)_theta$(theta_int).csv", data, ',')
            end
        end
    end 

end

function anopp2_noise(num_blades, v, omega, radii, dradii, chord, fn, fc, stationary_observer, theta)
    num_radial = length(radii)
    t0 = 0.0
    rot_axis = [0.0, 0.0, 1.0]
    blade_axis = [0.0, 1.0, 0.0]
    x0 = ANOPP2.A2_RK[cos(theta), 0.0, sin(theta)].*100.0.*12.0.*0.0254  # 100 ft in meters
    y0_hub = [0.0, 0.0, 0.0]  # m
    v0_hub = v.*rot_axis
    num_src_times = 256
    num_obs_times = 2*num_src_times

    # Now let's see if I can write out the ANOPP2 data.
    atm_conf = ANOPP2.writerspy.write_atm_files(density=rho, temperature=c0^2/(gam*R), pressure=rho*c0^2/gam)
    atm_t = ANOPP2.a2_atm_create(atm_conf)
    adl_confs = ANOPP2.writerspy.write_adl_files(
        num_blades=num_blades, num_radial=num_radial,
        num_src_times=num_src_times, aauc=area_over_chord_squared, omega=omega,
        radii=radii, dradii=dradii, chord=chord, fn=fn, fc=fc)
    surf_ts = [ANOPP2.a2_geo_create(adl_conf) for adl_conf in adl_confs]
    fp_conf = ANOPP2.writerspy.write_flightpath_files(rot_axis=rot_axis, blade_axis=blade_axis, v0_hub=v0_hub, y0_hub=y0_hub)
    fp_t = ANOPP2.a2_fp_create(fp_conf, atm_t)
    if stationary_observer
        # This just says the observer frame of reference is centered at the
        # global origin, and is not moving.
        obs_conf = ANOPP2.writerspy.write_observer_files(x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0])
    else
        # This just says the observer frame of reference is centered at the
        # global origin, and is moving with the same constant velocity as the hub.
        obs_conf = ANOPP2.writerspy.write_observer_files(x=[0.0, 0.0, 0.0], v=v0_hub)
    end
    obs_t = ANOPP2.a2_obs_create(obs_conf)
    ANOPP2.a2_obs_new_node(obs_t, x0)
    nnodes = ANOPP2.a2_obs_number_of_nodes(obs_t)
    f1a_conf = ANOPP2.writerspy.write_f1a_files(num_src_times=num_src_times, num_recep_times=num_obs_times, omega=omega)
    f1a_input_ts = vcat(surf_ts, atm_t)
    f1a_t, results_ts = ANOPP2.a2_exec_create_functional_module(f1a_conf, f1a_input_ts, obs_t)
    ANOPP2.a2_exec_execute_functional_module(f1a_t, atm_t, fp_t)
    t_monopole = Vector{Float64}[]
    t_dipole = Vector{Float64}[]
    p_monopole = Vector{Float64}[]
    p_dipole = Vector{Float64}[]
    for (j, results_t) in enumerate(results_ts)
        time, apth, nltype, cs, ife = ANOPP2.a2_obs_get_apth(obs_t, results_t, 1, 0)
        name = rstrip(ANOPP2.a2_obs_get_result_name(obs_t, results_t))
        if endswith(name, "Monopole")
            push!(t_monopole, time)
            push!(p_monopole, apth)
        elseif endswith(name, "Dipole")
            push!(t_dipole, time)
            push!(p_dipole, apth)
        end
    end

    p_monopole_total = p_monopole[1] .+ p_monopole[2]
    p_dipole_total = p_dipole[1] .+ p_dipole[2]

    return t_monopole[1], p_monopole_total, t_dipole[1], p_dipole_total

end

end # module
