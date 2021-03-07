function CompactSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ)
    ρ0 = op.rho
    c0 = op.asound
    r = section.r
    precone = rotor.precone
    T = promote_type(promote_type(typeof(r), typeof(precone)), typeof(θ))
    Np = out.Np
    Tp = out.Tp

    Λ = area_per_chord2*section.chord^2

    # Thinking about the geometry here:
    # y^     .            .
    #  |    .             .
    #  |   .              .
    #  |  .               .
    #  |φ.                .
    #  |.--------->x      .

    y0dot = @SVector [r*sin(precone), r*cos(precone)*cos(θ), r*cos(precone)*sin(θ)]
    y1dot = @SVector zeros(T, 3)
    y2dot = @SVector zeros(T, 3)
    y3dot = @SVector zeros(T, 3)

    T = promote_type(promote_type(typeof(Np), typeof(Tp)), typeof(θ))
    # The sign convention is a little tricky. We want the load on the fluid in
    # the blade coordinate system, which is a coordinate system that is rotating
    # and translating with blades, and assumes that the axis of rotation is at
    # the origin, aligned with the positive x axis. A positive rotation is
    # righthanded (so, say, if the blade is aligned with the y axis and rotating
    # with a positive rate, it is rotating toward the z axis) For the normal
    # loading, CCBlade gives a positive value when the load *on the blade* is in
    # the same direction as the axis of rotation (or opposite the freestream
    # velocity in the normal case). But we need the loading *on the fluid*,
    # which is in the opposite direction, hence we need to switch the sign on
    # Np.

    # For the circumferential loading, CCBlade gives a positive value when it
    # opposes the motion of the blade. So, in our hypothetical coordinate
    # system, let's say the blade is pointed along the y axis and rotating with
    # a positive rate. Then the blade is moving toward the positive z axis, and
    # since the circumferential loading opposes the blade motion, the load on
    # the blade would be pointed in the negative z direction. So that means the
    # load on the fluid would be in the positive z direction, so we don't need
    # to switch the sign.
    # 
    # Also decided to do the precone stuff. Need to check it somehow.
    f0dot = @SVector [-Np*cos(precone), -Tp*sin(θ) + Np*sin(precone)*cos(θ), Tp*cos(θ) + Np*sin(precone)*sin(θ)]
    f1dot = @SVector zeros(T, 3)

    return CompactSourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ)

end

function source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, num_blade_passes=2, steps_per_blade_pass=16)
    # Get the freestream parameters for each blade.
    Vinf = get_ccblade_Vinf(rotor, sections, ops)
    omega = get_ccblade_omega(rotor, sections, ops)

    # Need to know the radial spacing. (CCBlade doesn't use this—when
    # integrating stuff it uses the trapezoidal rule and passes in the radial
    # locations, and assumes that integrands go to zero at the hub and tip.)
    dradii = get_ccblade_dradii(rotor, sections)

    # Now we can get the transformations. Assume the rotor is traveling in the
    # positive x direction, with the first blade aligned with the
    # positive y axis. Rotor hub is at the origin.
    rot_axis = @SVector [1.0, 0.0, 0.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    v0_hub = Vinf.*rot_axis
    t0 = 0.0
    rot_trans = SteadyRotXTransformation(t0, omega, 0.0)
    const_vel_trans = ConstantVelocityTransformation(t0, y0_hub, v0_hub)

    # This is just an array of the angular offsets of each blade. First blade is
    # aligned with the y axis, next one is offset 2*pi/B radians, etc..
    B = rotor.B
    θs = 2*pi/B.*(0:(B-1))

    # Blade passing period (amount of time for one blade pass).
    bpp = 2*pi/omega/rotor.B

    # Get the time of each time step.
    src_time_range = num_blade_passes*bpp
    num_src_times = num_blade_passes*steps_per_blade_pass
    dt = src_time_range/(num_src_times - 1)
    src_times = t0 .+ (0:num_src_times-1).*dt

    # Get all the transformations.
    trans = compose.(src_times, Ref(const_vel_trans), Ref(rot_trans))

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    θs = reshape(θs, 1, 1, :)
    sections = reshape(sections, 1, :, 1)
    ops = reshape(ops, 1, :, 1)
    outputs = reshape(outputs, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Transform the source elements.
    ses = CompactSourceElement.(Ref(rotor), sections, ops, outputs, θs, dradii, area_per_chord2, src_times) .|> trans

    return ses
end

function get_ccblade_omega(rotor, sections, ops)
    radii = SingleFieldStructArray(sections, :r)
    Vys = SingleFieldStructArray(ops, :Vy)
    cos_precone = cos(rotor.precone)

    # Check: get the RPM from the radii and Vy.
    omegas = Vys./(radii*cos_precone)
    omega = first(omegas)
    # https://stackoverflow.com/a/51966131
    @assert all(x->x≈omega, omegas)

    return omega
end

function get_ccblade_Vinf(rotor, sections, ops)
    Vxs = SingleFieldStructArray(ops, :Vx)
    # Check: is the axial velocity is the same for each section?
    Vx = first(Vxs)
    # https://stackoverflow.com/a/51966131
    @assert all(x->x≈Vx, Vxs)
    # Should I worry about precone? I guess it wouldn't hurt.
    Vinf = Vx/cos(rotor.precone)

    return Vinf
end

function get_ccblade_dradii(rotor, sections)
    radii = SingleFieldStructArray(sections, :r)
    dradii = get_dradii(radii, rotor.Rhub, rotor.Rtip)
    return dradii
end
