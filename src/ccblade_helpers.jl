function _standard_ccblade_transform(rotor, sections, ops, period, num_src_times, positive_x_rotation)
    # Assume the rotor is traveling in the positive x direction, with the first
    # blade aligned with the positive y axis. Rotor hub is initially at the origin.
    rot_axis = @SVector [1.0, 0.0, 0.0]
    blade_axis = @SVector [0.0, 1.0, 0.0]
    y0_hub = @SVector [0.0, 0.0, 0.0]  # m
    t0 = 0.0

    # Get the time of each time step.
    dt = period/num_src_times
    src_times = t0 .+ (0:num_src_times-1).*dt

    # Get transformations for each blade element.
    cos_precone = cos(rotor.precone)
    r = SingleFieldStructArray(sections, Val{:r})
    Vx = SingleFieldStructArray(ops, Val{:Vx})
    Vy = SingleFieldStructArray(ops, Val{:Vy})
    if positive_x_rotation
        rot_trans = SteadyRotXTransformation.(t0, Vy./(r.*cos_precone), 0.0)  # size (num_radial,)
    else
        rot_trans = SteadyRotXTransformation.(t0, -Vy./(r.*cos_precone), 0.0)  # size (num_radial,)
    end
    const_vel_trans = ConstantVelocityTransformation.(t0, Ref(y0_hub), Ref(rot_axis).*Vx./cos_precone)  # size (num_radial,)

    # Reshape things to get broadcasting to work.
    rot_trans = reshape(rot_trans, 1, :)
    const_vel_trans = reshape(const_vel_trans, 1, :)

    # Now get all the transformations.
    trans = compose.(src_times, const_vel_trans, rot_trans)  # size (num_times, num_radial)

    return src_times, dt, trans
end

"""
    CompactSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ, positive_x_rotation=true)

Construct a source element to be used with the compact form of Farassat's formulation 1A from CCBlade objects.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.precone.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element.
- `area_per_chord2`: cross-sectional area divided by the chord squared of the element.
- `τ`: source time of the element.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function CompactSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ, positive_x_rotation)
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

    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = @SVector [r*spc, r*cpc*cθ, r*cpc*sθ]
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

    # So after all that, the takeaway is that we'll start out with a loading
    # vector [-Np, 0, Tp], then rotate it about the positive z-axis by an amount
    # `-precone`, then rotate it about the x axis by an amount `θ`.

    # But, what if I decide the blade is rotating around the negative x-axis?
    # The normal loading will still be in the same direction.
    # The precone and theta stuff can still work the same way.
    # I think the only thing that would switch is the circumferential loading.
    fn = -Np*cpc
    fr = Np*spc
    if positive_x_rotation
        fc = Tp
    else
        fc = -Tp
    end

    f0dot = @SVector [fn, cθ*fr - sθ*fc, sθ*fr + cθ*fc]
    T = eltype(f0dot)
    f1dot = @SVector zeros(T, 3)

    u = @SVector [spc, cpc*cθ, cpc*sθ]

    return CompactSourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ, u)

end

"""
    source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, area_per_chord2::Vector{AbstractFloat}, period, num_src_times, positive_x_rotation=true)

Construct and return an array of CompactSourceElement objects from CCBlade structs.

# Arguments
- `rotor`: CCBlade rotor object.
- `sections`: `Vector` of CCBlade section object.
- `ops`: `Vector` of CCBlade operating point.
- `outputs`::`Vector` of CCBlade output objects.
- `area_per_chord2`: cross-sectional area divided by the chord squared of the element at each CCBlade.section. Should be a Vector{AbstractFloat}, same length as `sections`, `ops`, `outputs`.
- `period`: length of the source time over which the returned source elements will evaluated.
- `num_src_times`: number of source times.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, period, num_src_times, positive_x_rotation)
    # Need to know the radial spacing. (CCBlade doesn't use this—when
    # integrating stuff [loading to get torque and thrust] it uses the
    # trapezoidal rule and passes in the radial locations, and assumes that
    # integrands go to zero at the hub and tip.) Kind of lame that I have to
    # calcluate it here, but whatever. Maybe I should use StaticArrays for this?
    # Ah, no, I don't know the length at compile time.
    dradii = get_ccblade_dradii(rotor, sections)

    # Get the transformation that will put the source elements in the "standard" CCBlade.jl reference frame (moving axially in the positive x axis direction, rotating about the positive x axis or negative x axis, first blade initially aligned with the positive y axis).
    src_times, dt, trans = _standard_ccblade_transform(rotor, sections, ops, period, num_src_times, positive_x_rotation)
    
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
    area_per_chord2 = reshape(area_per_chord2, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = CompactSourceElement.(Ref(rotor), sections, ops, outputs, θs, dradii, area_per_chord2, src_times, positive_x_rotation) .|> trans

    return ses
end


"""
    TBLTESourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)

Construct a source element to be used with the compact form of Farassat's formulation 1A from CCBlade objects.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.precone.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element, in meters.
- `τ`: source time of the element, in seconds.
- `Δτ`: source time duration, in seconds.
- `bl`: `AcousticAnalogies.AbstractBoundaryLayer`, needed for boundary layer properties.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function TBLTESourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)
    ρ0 = op.rho
    c0 = op.asound
    r = section.r
    precone = rotor.precone

    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    stwist, ctwist = sincos(section.theta + op.pitch)

    # The way this will work:
    #
    #   * We're going to assume that the rotor blade is moving axially in the positive x direction, and rotating about the x axis.
    #   * The "first" blade is initially aligned with the positive y axis.
    #   * Then we'll rotate the element about the positive z-axis by an amount `-precone` to acount for the precone.
    #   * Then we'll rotate the element about the positive x-axis by an amount `θ`.
    #
    # Those transformations will be applied to the blade element position, but also other vectors associated with the element.
    # In matrix form the precone transformation would be
    #
    # [ cos(-precone), -sin(-precone), 0 ]
    # [ sin(-precone),  cos(-precone), 0 ]
    # [      0,              0,        1 ]
    #
    # or equivalently
    #
    # [ cos(precone), sin(precone), 0 ]
    # [-sin(precone), cos(precone), 0 ]
    # [     0,             0,       1 ]
    #
    # And the θ rotation about the x axis would be
    #
    # [ 1,    0  ,    0    ]
    # [ 0, cos(θ), -sin(θ) ]
    # [ 0, sin(θ),  cos(θ) ]
    #
    # So, to get the position, we start out with a vector
    #
    # [0]
    # [r]
    # [0]
    #
    # And then multiply that by the precone rotation matrix and θ rotation matrix.
    #
    # [ cos(precone), sin(precone), 0 ] [0]   [ r*sin(precone) ]
    # [-sin(precone), cos(precone), 0 ] [r] = [ r*cos(precone) ]
    # [     0,             0,       1 ] [0]   [      0         ]
    #
    # [ 1,    0  ,    0    ] [ r*sin(precone) ]   [ r*sin(precone)        ]
    # [ 0, cos(θ), -sin(θ) ] [ r*cos(precone) ] = [ r*cos(precone)*cos(θ) ]
    # [ 0, sin(θ),  cos(θ) ] [      0         ]   [ r*cos(precone)*sin(θ) ]
    r = section.r
    y0dot = @SVector [r*spc, r*cpc*cθ, r*cpc*sθ]

    # In the blade-fixed frame, the source isn't moving, since the blade-fixed reference frame is moving with the source.
    T = eltype(y0dot)
    y1dot = @SVector zeros(T, 3)

    # Vx = op.Vx
    # u = out.u
    # Vy = op.Vy
    # v = out.v
    sphi, cphi = sincos(out.phi)
    Vx_plus_u = out.W*sphi
    Vy_minus_v = out.W*cphi

    # The `span_uvec` is a unit vector pointing from the hub to the tip, along the blade element's radial length.
    # So that's just the same as the position vector, but without the r factor.
    span_uvec = @SVector [spc, cpc*cθ, cpc*sθ]

    if positive_x_rotation
        # Now, what is the velocity of the fluid in the blade-fixed frame?
        # In our coordinate system, the rotor is rotating about the x axis, moving in the x axis direction.
        # So that means it appears that the axial freestream velocity Vx is in the negative x axis direction.
        # And the induced velocity `u` has the same sign convention as Vx.
        # So the total axial velocity of the fluid is `(-Vx - u)`.
        #
        # For the tangential velocity, we're imagining the blade is initially aligned with the y axis, rotating about the positive x axis.
        # So that means the blade is moving toward the z axis.
        # So, from the perspective of the blade, the `Vy` velocity is in the negative z axis direction.
        # But the sign convention for the induced tangential velocity `v` is that
        # it's positive when it opposes `Vy`, so the  total tangential velocity of the fluid is `(-Vy + v)`.
        #
        # So, finally, the fluid velocity vector we need to rotate is
        #
        # [-Vx - u ]
        # [    0   ]
        # [-Vy + v ]
        #
        # and when I do all that I get
        #
        # [ (-Vx - u)*cos(precone)                          ]
        # [ (Vx + u)*sin(precone)*cos(θ) - (-Vy + v)*sin(θ) ]
        # [ (Vx + u)*sin(precone)*sin(θ) + (-Vy + v)*cos(θ) ]
        # y1dot_fluid = @SVector [(-Vx - u)*cpc, (Vx + u)*spc*cθ - (-Vy + v)*sθ, (Vx + u)*spc*sθ + (-Vy + v)*cθ]
        y1dot_fluid = @SVector [-Vx_plus_u*cpc, Vx_plus_u*spc*cθ - (-Vy_minus_v)*sθ, Vx_plus_u*spc*sθ + (-Vy_minus_v)*cθ]

        # Finally the `chord_uvec` is a unit vector pointing from the leading edge to the trailing edge.
        # In our initial coordinate system (i.e., not accounting for the precone or
        # θ rotations) we're imagining the blade is rotating about the x axis,
        # aligned with the y axis, and so is moving in the direction of the z axis.
        # So that means if the twist is zero, then the leading edge is headed in the
        # z axis direction, and a vector pointing from leading edge to trailing edge
        # would be in the negative z axis direction. Then, to account for the twist,
        # we would rotate it about the positive y axis. And then do the usual
        # precone and θ rotations.
        # So, a rotation about the y axis is 
        #
        # [ cos(twist) 0 sin(twist) ]
        # [ 0          1       0    ]
        # [-sin(twist) 0 cos(twist) ]
        #
        # So, start with
        #
        # [ 0 ]
        # [ 0 ]
        # [-1 ]
        #
        # then
        #
        # [ cos(twist) 0 sin(twist) ] [ 0 ]   [-sin(twist) ]
        # [ 0          1       0    ] [ 0 ] = [    0       ]
        # [-sin(twist) 0 cos(twist) ] [-1 ]   [-cos(twist) ]
        #
        # Now do the precone transformation
        #
        # [ cos(precone), sin(precone), 0 ] [-sin(twist) ]   [-sin(twist)*cos(precone) ]
        # [-sin(precone), cos(precone), 0 ] [    0       ] = [ sin(twist)*sin(precone) ]
        # [     0,             0,       1 ] [-cos(twist) ]   [-cos(twist)              ]
        #
        # Finally do the θ transformation
        #
        # [ 1,    0  ,    0    ] [-sin(twist)*cos(precone) ]   [-sin(twist)*cos(precone)                            ]
        # [ 0, cos(θ), -sin(θ) ] [ sin(twist)*sin(precone) ] = [ sin(twist)*sin(precone)*cos(θ) + cos(twist)*sin(θ) ]
        # [ 0, sin(θ),  cos(θ) ] [-cos(twist)              ]   [ sin(twist)*sin(precone)*sin(θ) - cos(twist)*cos(θ) ]
        chord_uvec = @SVector [-stwist*cpc, stwist*spc*cθ + ctwist*sθ, stwist*spc*sθ - ctwist*cθ]

    else

        # But, what if I want to assume that the blade is rotating in the opposite direction, i.e., about the negative x axis?
        # For the velocity, the direction of the axial velocity is unchanged: we're still moving in the positive x direction, so the axial velocity from the perspective of the blade element will be in the negative x direction.
        # So the total axial velocity of the fluid is `(-Vx - u)`.
        #
        # For the tangential velocity, we're rotating about the negative x axis now, so since the blade is initially aligned with the y axis, it is moving in the negative z direction.
        # So that means the freestream tangential velocity appears to be in the opposite direction, aka the positive z direction.
        # But the induced tangential velocity is in the opposite direction of the freestream tangential velocity, so the total velocity in the tangential direction is `(Vy - v)`.
        # 
        # So the fluid velocity vector we want to rotate is
        #
        # [-Vx - u ]
        # [    0   ]
        # [ Vy - v ]
        #
        # The theta and precone stuff doesn't change, so we'll do all the same stuff.
        # First we do the precone:
        #
        # [ cos(precone), sin(precone), 0 ] [-Vx - u]   [ (-Vx - u)*cos(precone)    ]   [ (-Vx - u)*cos(precone) ]
        # [-sin(precone), cos(precone), 0 ] [    0  ] = [ (-Vx - u)*(-sin(precone)) ] = [ ( Vx + u)*sin(precone)  ]
        # [     0,             0,       1 ] [ Vy - v]   [      Vy - v               ]   [   Vy - v                ]
        #
        # then do the theta rotation.
        #
        # [ 1,    0  ,    0    ] [ (-Vx - u)*cos(precone) ]   [ (-Vx - u)*cos(precone)                          ]
        # [ 0, cos(θ), -sin(θ) ] [ ( Vx + u)*sin(precone) ] = [ ( Vx + u)*sin(precone)*cos(θ) - (Vy - v)*sin(θ) ]
        # [ 0, sin(θ),  cos(θ) ] [   Vy - v               ]   [ ( Vx + u)*sin(precone)*sin(θ) + (Vy - v)*cos(θ) ]
        # y1dot_fluid = @SVector [(-Vx - u)*cpc, (Vx + u)*spc*cθ - (Vy - v)*sθ, (Vx + u)*spc*sθ + (Vy - v)*cθ]
        y1dot_fluid = @SVector [-Vx_plus_u*cpc, Vx_plus_u*spc*cθ - Vy_minus_v*sθ, Vx_plus_u*spc*sθ + Vy_minus_v*cθ]
        #
        # That should be the same thing as the opposite case, but with the sign on (Vy - v) switched.
        # Yep, good.
        #
        # For the chord_uvec, I want to start with a unit vector pointing in the positive z axis, then do a negative-twist rotation about the positive y axis.
        # So start with
        #
        # [ 0 ]
        # [ 0 ]
        # [ 1 ]
        #
        # then
        #
        # [ cos(-twist) 0 sin(-twist) ] [ 0 ]   [ sin(-twist) ]
        # [ 0            1       0    ] [ 0 ] = [      0      ]
        # [-sin(-twist) 0 cos(-twist) ] [ 1 ]   [ cos(-twist) ]
        #
        # Now do the precone transformation
        #
        # [ cos(precone), sin(precone), 0 ] [ sin(-twist) ]   [ sin(-twist)*cos(precone) ]
        # [-sin(precone), cos(precone), 0 ] [     0       ] = [-sin(-twist)*sin(precone) ]
        # [     0,             0,       1 ] [ cos(-twist) ]   [ cos(-twist)              ]
        #
        # Finally do the θ transformation
        #
        # [ 1,    0  ,    0    ] [ sin(-twist)*cos(precone) ]   [ sin(-twist)*cos(precone)                             ]
        # [ 0, cos(θ), -sin(θ) ] [-sin(-twist)*sin(precone) ] = [-sin(-twist)*sin(precone)*cos(θ) - cos(-twist)*sin(θ) ]
        # [ 0, sin(θ),  cos(θ) ] [ cos(-twist)              ]   [-sin(-twist)*sin(precone)*sin(θ) + cos(-twist)*cos(θ) ]
        #
        # Now handle the `-twist`,
        #
        # [ sin(-twist)*cos(precone)                             ]   [-sin(twist)*cos(precone)                            ]
        # [-sin(-twist)*sin(precone)*cos(θ) - cos(-twist)*sin(θ) ] = [ sin(twist)*sin(precone)*cos(θ) - cos(twist)*sin(θ) ]
        # [-sin(-twist)*sin(precone)*sin(θ) + cos(-twist)*cos(θ) ]   [ sin(twist)*sin(precone)*sin(θ) + cos(twist)*cos(θ) ]
        
        chord_uvec = @SVector [-stwist*cpc, stwist*spc*cθ - ctwist*sθ, stwist*spc*sθ + ctwist*cθ]
    end

    c0 = op.asound
    nu = op.mu/op.rho
    chord = section.chord

    chord_cross_span_to_get_top_uvec = positive_x_rotation
    return TBLTESourceElement(c0, nu, Δr, chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    tblte_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, bls::Vector{AbstractBoundaryLayer}, period, num_src_times, positive_x_rotation)

Construct and return an array of CompactSourceElement objects from CCBlade structs.

# Arguments
- `rotor`: CCBlade rotor object.
- `sections`: `Vector` of CCBlade section object.
- `ops`: `Vector` of CCBlade operating point.
- `outputs`: `Vector` of CCBlade output objects.
- `bls`::`Vector` of boundary layer `AbstractBoundaryLayer` `structs`.
- `period`: length of the source time over which the returned source elements will evaluated.
- `num_src_times`: number of source times.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function tblte_source_elements_ccblade(rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
    # Need to know the radial spacing. (CCBlade doesn't use this—when
    # integrating stuff [loading to get torque and thrust] it uses the
    # trapezoidal rule and passes in the radial locations, and assumes that
    # integrands go to zero at the hub and tip.) Kind of lame that I have to
    # calcluate it here, but whatever. Maybe I should use StaticArrays for this?
    # Ah, no, I don't know the length at compile time.
    dradii = get_ccblade_dradii(rotor, sections)

    # Get the transformation that will put the source elements in the "standard" CCBlade.jl reference frame (moving axially in the positive x axis direction, rotating about the positive x axis, first blade initially aligned with the positive y axis).
    src_times, dt, trans = _standard_ccblade_transform(rotor, sections, ops, period, num_src_times, positive_x_rotation)

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
    bls = reshape(bls, 1, :, 1)
    src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = TBLTESourceElement.(Ref(rotor), sections, ops, outputs, θs, dradii, src_times, Ref(dt), bls, positive_x_rotation) .|> trans

    return ses
end

"""
    get_ccblade_dradii(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section})

Construct and return a Vector of the lengths of each CCBlade section.
"""
function get_ccblade_dradii(rotor, sections)
    radii = SingleFieldStructArray(sections, Val{:r})
    dradii = get_dradii(radii, rotor.Rhub, rotor.Rtip)
    return dradii
end
