"""
    CompactSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ)

Construct a source element to be used with the compact form of Farassat's formulation 1A from CCBlade objects.

# Arguments
- 
"""
function CompactSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ)
    ρ0 = op.rho
    c0 = op.asound
    r = section.r
    precone = rotor.precone
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
    T = eltype(y0dot)
    y1dot = @SVector zeros(T, 3)
    y2dot = @SVector zeros(T, 3)
    y3dot = @SVector zeros(T, 3)

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
    # load on the fluid would be in the positive z direction, and we don't need
    # to switch the sign.
    fn = -Np*cos(precone)
    fr = Np*sin(precone)
    fc = Tp

    f0dot = @SVector [fn, cos(θ)*fr - sin(θ)*fc, sin(θ)*fr + cos(θ)*fc]
    T = eltype(f0dot)
    f1dot = @SVector zeros(T, 3)

    return CompactSourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ)

end

function source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, period, num_src_times)
    # Need to know the radial spacing. (CCBlade doesn't use this—when
    # integrating stuff (torque and thrust) it uses the trapezoidal rule and
    # passes in the radial locations, and assumes that integrands go to zero at
    # the hub and tip.) Kind of lame that I have to calcluate it here, but
    # whatever.
    dradii = get_ccblade_dradii(rotor, sections)

    # Assume the rotor is traveling in the positive x direction, with the first
    # blade aligned with the positive y axis. Rotor hub is at the origin.
    rot_axis = @SVector [1.0, 0.0, 0.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    t0 = 0.0

    # Get the time of each time step.
    dt = period/(num_src_times - 1)
    src_times = t0 .+ (0:num_src_times-1).*dt

    # Get transformations for each blade element.
    cos_precone = cos(rotor.precone)
    r = SingleFieldStructArray(sections, :r)
    Vx = SingleFieldStructArray(ops, :Vx)
    Vy = SingleFieldStructArray(ops, :Vy)
    rot_trans = SteadyRotXTransformation.(t0, Vy./(r.*cos_precone), 0.0)  # size (num_radial,)
    const_vel_trans = ConstantVelocityTransformation.(t0, Ref(y0_hub), Ref(rot_axis).*Vx./cos_precone)  # size (num_radial,)

    # Reshape things to get broadcasting to work.
    rot_trans = reshape(rot_trans, 1, :)
    const_vel_trans = reshape(const_vel_trans, 1, :)

    # Now get all the transformations.
    trans = compose.(src_times, const_vel_trans, rot_trans)  # size (num_times, num_radial)
    
    # This is just an array of the angular offsets of each blade. First blade is
    # aligned with the y axis, next one is offset 2*pi/B radians, etc..
    num_blades = rotor.B
    θs = 2*pi/num_blades.*(0:(num_blades-1))

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    trans = reshape(trans, size(trans)..., 1)
    θs = reshape(θs, 1, 1, :)
    sections = reshape(sections, 1, :, 1)
    ops = reshape(ops, 1, :, 1)
    outputs = reshape(outputs, 1, :, 1)
    dradii = reshape(dradii, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = CompactSourceElement.(Ref(rotor), sections, ops, outputs, θs, dradii, area_per_chord2, src_times) .|> trans

    return ses
end

function get_ccblade_dradii(rotor, sections)
    radii = SingleFieldStructArray(sections, :r)
    dradii = get_dradii(radii, rotor.Rhub, rotor.Rtip)
    return dradii
end
