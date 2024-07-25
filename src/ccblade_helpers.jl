import FillArrays: getindex_value

# Normal implementation is this from FillArrays.jl:
#
# @inline getindex_value(F::Fill) = F.value
#
# But that breaks with CCBlade.
# Why?
# Because, since a FillArrays.Fill <: AbstractArray, it calls this:
#
# function Base.getproperty(obj::AbstractVector{<:Section}, sym::Symbol)
#     return getfield.(obj, sym)
# end
#
# which eventually calls FillArrays.getindex_value again, leading to recursion and a stack overflow.
#
# This is type piracy :-(.
# But it may also be type piracy to extend Base.getproperty in CCBlade.jl, since CCBlade.jl doesn't own Base.getproperty or AbstractVector.
@inline getindex_value(F::Fill{<:Union{CCBlade.Section,CCBlade.OperatingPoint,CCBlade.Outputs}}) = getfield(F, :value)

# Normal implementation is this from FillArrays.jl:
#
# @inline axes(F::Fill) = F.axes
#
# But that hits 
#
# function Base.getproperty(obj::AbstractVector{<:Section}, sym::Symbol)
#     return getfield.(obj, sym)
# end
#
# from CCBlade.
# This is type piracy :-(.
# But it may also be type piracy to extend Base.getproperty in CCBlade.jl, since CCBlade.jl doesn't own Base.getproperty or AbstractVector.
@inline Base.axes(F::Fill{<:Union{CCBlade.Section,CCBlade.OperatingPoint,CCBlade.Outputs}}) = getfield(F, :axes)

function _standard_ccblade_transform(rotor::CCBlade.Rotor, sections::AbstractVector{<:CCBlade.Section}, ops::AbstractVector{<:CCBlade.OperatingPoint}, period, num_src_times, positive_x_rotation)
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
    # r = SingleFieldStructArray(sections, Val{:r})
    # Vx = SingleFieldStructArray(ops, Val{:Vx})
    # Vy = SingleFieldStructArray(ops, Val{:Vy})
    r = mapview(:r, sections)
    Vx = mapview(:Vx, ops)
    Vy = mapview(:Vy, ops)
    if positive_x_rotation
        rot_trans = SteadyRotXTransformation.(t0, Vy./(r.*cos_precone), 0.0)  # size (num_radial,)
    else
        rot_trans = SteadyRotXTransformation.(t0, -Vy./(r.*cos_precone), 0.0)  # size (num_radial,)
    end
    const_vel_trans = ConstantVelocityTransformation.(t0, Ref(y0_hub), Ref(rot_axis).*Vx./cos_precone)  # size (num_radial,)

    # Reshape things to get broadcasting to work.
    rot_trans_rs = reshape(rot_trans, 1, :)
    const_vel_trans_rs = reshape(const_vel_trans, 1, :)

    # Now get all the transformations.
    trans = compose.(src_times, const_vel_trans_rs, rot_trans_rs)  # size (num_times, num_radial)

    return src_times, dt, trans
end

"""
    CompactF1ASourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ, positive_x_rotation=true)

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
function CompactF1ASourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, area_per_chord2, τ, positive_x_rotation)
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

    return CompactF1ASourceElement(ρ0, c0, Δr, Λ, y0dot, y1dot, y2dot, y3dot, f0dot, f1dot, τ, u)

end

"""
    f1a_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, area_per_chord2::Vector{AbstractFloat}, period, num_src_times, positive_x_rotation)

Construct and return an array of CompactF1ASourceElement objects from CCBlade structs.

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
function f1a_source_elements_ccblade(rotor, sections, ops, outputs, area_per_chord2, period, num_src_times, positive_x_rotation)
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
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    # trans_rs = reshape(trans, size(trans)..., 1)
    θs_rs = reshape(θs, 1, 1, :)
    sections_rs = reshape(sections, 1, :, 1)
    ops_rs = reshape(ops, 1, :, 1)
    outputs_rs = reshape(outputs, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    area_per_chord2_rs = reshape(area_per_chord2, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = CompactF1ASourceElement.(Ref(rotor), sections_rs, ops_rs, outputs_rs, θs_rs, dradii_rs, area_per_chord2_rs, src_times, positive_x_rotation) .|> trans

    return ses
end

function _get_position_velocity_span_uvec_chord_uvec(theta, precone, pitch, r, θ, W, phi, positive_x_rotation)
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    stwist, ctwist = sincos(theta + pitch)

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
    y0dot = @SVector [r*spc, r*cpc*cθ, r*cpc*sθ]

    # In the blade-fixed frame, the source isn't moving, since the blade-fixed reference frame is moving with the source.
    y1dot = @SVector zeros(eltype(y0dot), 3)

    # Vx = op.Vx
    # u = out.u
    # Vy = op.Vy
    # v = out.v
    sphi, cphi = sincos(phi)
    Vx_plus_u = W*sphi
    Vy_minus_v = W*cphi

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

    chord_cross_span_to_get_top_uvec = positive_x_rotation
    return y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec
end


"""
    TBLTESourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)

Construct a source element to be used to predict turbulent boundary layer-trailing edge (TBLTE) noise.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
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
    return TBLTESourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(rotor, section, op, out, θ, Δr, τ, Δτ, bl, positive_x_rotation)
end

function TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation) where {TDirect,TUInduction,TMachCorrection,TDoppler}
   
    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    nu = op.mu/op.rho

    return TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(op.asound, nu, Δr, section.chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    tblte_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, bls::Vector{AbstractBoundaryLayer}, period, num_src_times, positive_x_rotation)

Construct and return an array of TBLTESourceElement objects from CCBlade structs.

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
    return tblte_source_elements_ccblade(BrooksBurleyDirectivity, true, PrandtlGlauertMachCorrection, true, rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
end

function tblte_source_elements_ccblade(TDirect::Type{<:AbstractDirectivity}, TUInduction::Bool, TMachCorrection::Type{<:AbstractMachCorrection}, TDoppler::Bool, rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
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
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    # trans_rs = reshape(trans, size(trans)..., 1)
    θs_rs = reshape(θs, 1, 1, :)
    sections_rs = reshape(sections, 1, :, 1)
    ops_rs = reshape(ops, 1, :, 1)
    outputs_rs = reshape(outputs, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    bls_rs = reshape(bls, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = TBLTESourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}.(Ref(rotor), sections_rs, ops_rs, outputs_rs, θs_rs, dradii_rs, src_times, Ref(dt), bls_rs, positive_x_rotation) .|> trans

    return ses
end

"""
    LBLVSSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)

Construct a source element to be used to predict laminary boundary layer-vortex shedding (LBLVS) noise.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
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
function LBLVSSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)
    return LBLVSSourceElement{BrooksBurleyDirectivity,true,true}(rotor, section, op, out, θ, Δr, τ, Δτ, bl, positive_x_rotation)
end

function LBLVSSourceElement{TDirect,TUInduction,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation) where {TDirect,TUInduction,TDoppler}

    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    nu = op.mu/op.rho

    return LBLVSSourceElement{TDirect,TUInduction,TDoppler}(op.asound, nu, Δr, section.chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    lblvs_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, bls::Vector{AbstractBoundaryLayer}, period, num_src_times, positive_x_rotation)

Construct and return an array of LBLVSSourceElement objects from CCBlade structs.

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
function lblvs_source_elements_ccblade(rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
    return lblvs_source_elements_ccblade(BrooksBurleyDirectivity, true, true, rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
end

function lblvs_source_elements_ccblade(TDirect::Type{<:AbstractDirectivity}, TUInduction::Bool, TDoppler::Bool, rotor, sections, ops, outputs, bls, period, num_src_times, positive_x_rotation)
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
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    # trans_rs = reshape(trans, size(trans)..., 1)
    θs_rs = reshape(θs, 1, 1, :)
    sections_rs = reshape(sections, 1, :, 1)
    ops_rs = reshape(ops, 1, :, 1)
    outputs_rs = reshape(outputs, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    bls_rs = reshape(bls, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = LBLVSSourceElement{TDirect,TUInduction,TDoppler}.(Ref(rotor), sections_rs, ops_rs, outputs_rs, θs_rs, dradii_rs, src_times, Ref(dt), bls_rs, positive_x_rotation) .|> trans

    return ses
end

"""
    TipVortexSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation)

Construct a source element to be used to predict tip vortex noise.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element, in meters.
- `τ`: source time of the element, in seconds.
- `Δτ`: source time duration, in seconds.
- `bl`: `AcousticAnalogies.AbstractBoundaryLayer`, needed for boundary layer properties.
- `blade_tip`: `AcousticAnalogies.AbstractBladeTip`
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function TipVortexSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation)
    return TipVortexSourceElement{BrooksBurleyDirectivity,true,true}(rotor, section, op, out, θ, Δr, τ, Δτ, bl, blade_tip, positive_x_rotation)
end

function TipVortexSourceElement{TDirect,TUInduction,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation) where {TDirect,TUInduction,TDoppler}

    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    return TipVortexSourceElement{TDirect,TUInduction,TDoppler}(op.asound, Δr, section.chord, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

"""
    tip_vortex_source_elements_ccblade(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, output::CCBlade.Outputs, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, period, num_src_times, positive_x_rotation)

Construct and return an array of TipVortexSourceElement objects from CCBlade structs.

Note that unlike the other `*_source_elements_ccblade` functions, `tip_vortex_source_elements_ccblade` expects scalar arguments instead of vectors for `section`, `op`, etc. as a blade only has one tip.

# Arguments
- `rotor`: CCBlade rotor object.
- `section`: CCBlade section object at the blade tip.
- `op`: CCBlade operating point object at the blade tip.
- `output`: CCBlade output object at the blade tip.
- `Δr`: radial spacing.
- `bl`:: Boundary layer `struct` at the blade tip.
- `blade_tip`: `AcousticAnalogies.AbstractBladeTip`
- `period`: length of the source time over which the returned source elements will evaluated.
- `num_src_times`: number of source times.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function tip_vortex_source_elements_ccblade(rotor, section, op, output, Δr, bl, blade_tip, period, num_src_times, positive_x_rotation)
    return tip_vortex_source_elements_ccblade(BrooksBurleyDirectivity, true, true, rotor, section, op, output, Δr, bl, blade_tip, period, num_src_times, positive_x_rotation)
end

function tip_vortex_source_elements_ccblade(TDirect::Type{<:AbstractDirectivity}, TUInduction::Bool, TDoppler::Bool, rotor, section, op, output, Δr, bl, blade_tip, period, num_src_times, positive_x_rotation)
    # Ugh, hate doing this.
    # Wish there was a way to make a allocation-free array-like thingy from a scaler.
    # But I doubt it makes any difference.
    # sections = [section]
    # ops = [op]
    # Good news!
    # Learned about the FillArrays.jl package.
    sections = Fill(section, 1)
    ops = Fill(op, 1)
    # But that breaks with CCBlade.jl.
    # So back to 1D arrays.
    # sections = [section]
    # ops = [op]

    # Get the transformation that will put the source elements in the "standard" CCBlade.jl reference frame (moving axially in the positive x axis direction, rotating about the positive x axis, first blade initially aligned with the positive y axis).
    src_times, dt, trans = _standard_ccblade_transform(rotor, sections, ops, period, num_src_times, positive_x_rotation)

    # This is just an array of the angular offsets of each blade. First blade is
    # aligned with the y axis, next one is offset 2*pi/B radians, etc..
    num_blades = rotor.B
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times, num_radial, num_blades).
    # But this will really be (num_times, 1, num_blades).
    # trans_rs = reshape(trans, size(trans)..., 1)
    θs_rs = reshape(θs, 1, 1, :)
    # sections_rs = reshape(sections, 1, 1)
    # ops_rs = reshape(ops, 1, 1)
    # outputs = reshape(outputs, 1, :, 1)
    # dradii = reshape(dradii, 1, :, 1)
    # bls = reshape(bls, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    # ses = TipVortexSourceElement.(Ref(rotor), sections_rs, ops_rs, Ref(output), θs_rs, Ref(Δr), src_times, Ref(dt), Ref(bl), positive_x_rotation) .|> trans_rs
    # So Θs_rs has size (1, 1, num_blades), src_times has size (num_src_times,), trans has size (num_src_times, 1).
    # So that should all work out.
    ses = TipVortexSourceElement{TDirect,TUInduction,TDoppler}.(Ref(rotor), Ref(section), Ref(op), Ref(output), θs_rs, Ref(Δr), src_times, Ref(dt), Ref(bl), Ref(blade_tip), positive_x_rotation) .|> trans

    return ses
end

"""
    TEBVSSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)

Construct a source element to be used to predict trailing edge bluntness-vortex shedding (TEBVS) noise.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element, in meters.
- `h`: trailing edge thickness (m)
- `Psi`: solid angle between the blade surfaces immediately upstream of the trailing edge (rad)
- `τ`: source time of the element, in seconds.
- `Δτ`: source time duration, in seconds.
- `bl`: `AcousticAnalogies.AbstractBoundaryLayer`, needed for boundary layer properties.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function TEBVSSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)
    return TEBVSSourceElement{BrooksBurleyDirectivity,true,true}(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation)
end

function TEBVSSourceElement{TDirect,TUInduction,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation) where {TDirect,TUInduction,TDoppler}

    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    nu = op.mu/op.rho

    return TEBVSSourceElement{TDirect,TUInduction,TDoppler}(op.asound, nu, Δr, section.chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    tebvs_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, hs, Psis, bls::Vector{AbstractBoundaryLayer}, period, num_src_times, positive_x_rotation)

Construct and return an array of TEBVSSourceElement objects from CCBlade structs.

# Arguments
- `rotor`: CCBlade rotor object.
- `sections`: `Vector` of CCBlade section object.
- `ops`: `Vector` of CCBlade operating point.
- `outputs`: `Vector` of CCBlade output objects.
- `hs`: `Vector` of trailing edge thicknesses
- `Psis`: `Vector` of solid angles between the blade surfaces immediately upstream of the trailing edge (rad)
- `bls`::`Vector` of boundary layer `AbstractBoundaryLayer` `structs`.
- `period`: length of the source time over which the returned source elements will evaluated.
- `num_src_times`: number of source times.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function tebvs_source_elements_ccblade(rotor, sections, ops, outputs, hs, Psis, bls, period, num_src_times, positive_x_rotation)
    return tebvs_source_elements_ccblade(BrooksBurleyDirectivity, true, true, rotor, sections, ops, outputs, hs, Psis, bls, period, num_src_times, positive_x_rotation)
end

function tebvs_source_elements_ccblade(TDirect::Type{<:AbstractDirectivity}, TUInduction::Bool, TDoppler::Bool, rotor, sections, ops, outputs, hs, Psis, bls, period, num_src_times, positive_x_rotation)
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
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times,
    # num_radial, num_blades).
    # trans_rs = reshape(trans, size(trans)..., 1)
    θs_rs = reshape(θs, 1, 1, :)
    sections_rs = reshape(sections, 1, :, 1)
    ops_rs = reshape(ops, 1, :, 1)
    outputs_rs = reshape(outputs, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    hs_rs = reshape(hs, 1, :, 1)
    Psis_rs = reshape(Psis, 1, :, 1)
    bls_rs = reshape(bls, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # Construct and transform the source elements.
    ses = TEBVSSourceElement{TDirect,TUInduction,TDoppler}.(Ref(rotor), sections_rs, ops_rs, outputs_rs, θs_rs, dradii_rs, hs_rs, Psis_rs, src_times, Ref(dt), bls_rs, positive_x_rotation) .|> trans

    return ses
end

"""
    CombinedNoTipBroadbandSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE), laminar boundary layer-vortex shedding (LBLVS) noise, and trailing edge bluntness-vortex shedding (TEBVS) noise using the BPM/Brooks and Burley method from CCBlade structs.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element, in meters.
- `h`: trailing edge thickness (m)
- `Psi`: solid angle between the blade surfaces immediately upstream of the trailing edge (rad)
- `τ`: source time of the element, in seconds.
- `Δτ`: source time duration, in seconds.
- `bl`: `AcousticAnalogies.AbstractBoundaryLayer`, needed for boundary layer properties.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function CombinedNoTipBroadbandSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation)
    return CombinedNoTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, positive_x_rotation)
end

function CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, positive_x_rotation) where {TDirect,TUInduction,TMachCorrection,TDoppler}

    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    nu = op.mu/op.rho

    return CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(op.asound, nu, Δr, section.chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, chord_cross_span_to_get_top_uvec)
end

"""
    CombinedWithTipBroadbandSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation)

Construct a source element for predicting turbulent boundary layer-trailing edge (TBLTE), laminar boundary layer-vortex shedding (LBLVS) noise, trailing edge bluntness-vortex shedding (TEBVS), and tip vortex noise using the BPM/Brooks and Burley method from CCBlade structs.

The source element's position is calculated from `section.r`, `rotor.precone`, and the `θ` argument using
```julia
    sθ, cθ = sincos(θ)
    spc, cpc = sincos(precone)
    y0dot = [r*spc, r*cpc*cθ, r*cpc*sθ]
```
where `y0dot` is the position of the source element.

# Arguments
- `rotor::CCBlade.Rotor`: CCBlade rotor object, needed for the precone angle.
- `section::CCBlade.Section`: CCBlade section object, needed for the radial location and chord length of the element.
- `op::CCBlade.OperatingPoint`: CCBlade operating point, needed for atmospheric properties.
- `out::CCBlade.Outputs`: CCBlade outputs object, needed for the loading.
- `θ`: polar coordinate of the element, in radians.
- `Δr`: length of the element, in meters.
- `h`: trailing edge thickness (m)
- `Psi`: solid angle between the blade surfaces immediately upstream of the trailing edge (rad)
- `τ`: source time of the element, in seconds.
- `Δτ`: source time duration, in seconds.
- `bl`: `AcousticAnalogies.AbstractBoundaryLayer`, needed for boundary layer properties.
- `blade_tip`: Blade tip struct, i.e. an AbstractBladeTip.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function CombinedWithTipBroadbandSourceElement(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation)
    return CombinedWithTipBroadbandSourceElement{BrooksBurleyDirectivity,true,PrandtlGlauertMachCorrection,true}(rotor, section, op, out, θ, Δr, h, Psi, τ, Δτ, bl, blade_tip, positive_x_rotation)
end

function CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(rotor::CCBlade.Rotor, section::CCBlade.Section, op::CCBlade.OperatingPoint, out::CCBlade.Outputs, θ, Δr, h, Psi, τ, Δτ, bl::AbstractBoundaryLayer, blade_tip::AbstractBladeTip, positive_x_rotation) where {TDirect,TUInduction,TMachCorrection,TDoppler}

    y0dot, y1dot, y1dot_fluid, span_uvec, chord_uvec, chord_cross_span_to_get_top_uvec = _get_position_velocity_span_uvec_chord_uvec(
        section.theta, rotor.precone, op.pitch, section.r, θ, out.W, out.phi, positive_x_rotation)

    nu = op.mu/op.rho

    return CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}(op.asound, nu, Δr, section.chord, h, Psi, y0dot, y1dot, y1dot_fluid, τ, Δτ, span_uvec, chord_uvec, bl, blade_tip, chord_cross_span_to_get_top_uvec)
end

"""
    combined_broadband_source_elements_ccblade(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section}, ops::Vector{CCBlade.OperatingPoint}, outputs::Vector{CCBlade.Outputs}, hs::Vector{Float64}, Psis::Vector{Float64}, bls::Vector{AbstractBoundaryLayer}, blade_tip::AbstractBladeTip, period, num_src_times, positive_x_rotation)

Construct and return an array of broadband prediction source element objects from CCBlade structs.

# Arguments
- `rotor`: CCBlade rotor object.
- `sections`: `Vector` of CCBlade section object.
- `ops`: `Vector` of CCBlade operating point.
- `outputs`: `Vector` of CCBlade output objects.
- `hs`: `Vector` of trailing edge thicknesses (m)
- `Psis`: `Vector` of solid angles between the blade surfaces immediately upstream of the trailing edge (rad)
- `bls`::`Vector` of boundary layer `AbstractBoundaryLayer` `structs`.
- `blade_tip`: Blade tip struct, i.e. an AbstractBladeTip.
- `period`: length of the source time over which the returned source elements will evaluated.
- `num_src_times`: number of source times.
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise.
"""
function combined_broadband_source_elements_ccblade(rotor, sections, ops, outputs, hs, Psis, bls, blade_tip, period, num_src_times, positive_x_rotation)
    return combined_broadband_source_elements_ccblade(BrooksBurleyDirectivity, true, PrandtlGlauertMachCorrection, true, rotor, sections, ops, outputs, hs, Psis, bls, blade_tip, period, num_src_times, positive_x_rotation)
end

function combined_broadband_source_elements_ccblade(TDirect::Type{<:AbstractDirectivity}, TUInduction::Bool, TMachCorrection::Type{<:AbstractMachCorrection}, TDoppler::Bool, rotor, sections, ops, outputs, hs, Psis, bls::AbstractVector{<:AbstractBoundaryLayer}, blade_tip, period, num_src_times, positive_x_rotation)
    # Need to know the radial spacing. (CCBlade doesn't use this—when
    # integrating stuff [loading to get torque and thrust] it uses the
    # trapezoidal rule and passes in the radial locations, and assumes that
    # integrands go to zero at the hub and tip.) Kind of lame that I have to
    # calcluate it here, but whatever. Maybe I should use StaticArrays for this?
    # Ah, no, I don't know the length at compile time.
    dradii = get_ccblade_dradii(rotor, sections)

    # Get the transformation that will put the source elements in the "standard" CCBlade.jl reference frame (moving axially in the positive x axis direction, rotating about the positive x axis, first blade initially aligned with the positive y axis).
    # Will be size (num_times, num_radial), so we'll need to adjust for the no tip/with tip stuff.
    src_times, dt, trans = _standard_ccblade_transform(rotor, sections, ops, period, num_src_times, positive_x_rotation)

    # This is just an array of the angular offsets of each blade. First blade is
    # aligned with the y axis, next one is offset 2*pi/B radians, etc..
    num_blades = rotor.B
    θs = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)

    # Reshape for broadcasting. Goal is to make everything work for a size of (num_times, num_radial, num_blades).
    θs_rs = reshape(θs, 1, 1, :)
    sections_rs = reshape(sections, 1, :, 1)
    ops_rs = reshape(ops, 1, :, 1)
    outputs_rs = reshape(outputs, 1, :, 1)
    dradii_rs = reshape(dradii, 1, :, 1)
    hs_rs = reshape(hs, 1, :, 1)
    Psis_rs = reshape(Psis, 1, :, 1)
    bls_rs = reshape(bls, 1, :, 1)
    # src_times = reshape(src_times, :, 1, 1)  # This one isn't necessary.

    # So, I want to create some structs for all the non-blade tip elements, and the blade tip elements.
    # So I just need to slice things appropriately.
    sections_rs_no_tip = @view sections_rs[:, begin:end-1, :]
    ops_rs_no_tip = @view ops_rs[:, begin:end-1, :]
    outputs_rs_no_tip = @view outputs_rs[:, begin:end-1, :]
    dradii_rs_no_tip = @view dradii_rs[:, begin:end-1, :]
    hs_rs_no_tip = @view hs_rs[:, begin:end-1, :]
    Psis_rs_no_tip = @view Psis_rs[:, begin:end-1, :]
    bls_rs_no_tip = @view bls_rs[:, begin:end-1, :]
    trans_no_tip = @view trans[:, begin:end-1]

    sections_rs_with_tip = @view sections_rs[:, end:end, :]
    ops_rs_with_tip = @view ops_rs[:, end:end, :]
    outputs_rs_with_tip = @view outputs_rs[:, end:end, :]
    dradii_rs_with_tip = @view dradii_rs[:, end:end, :]
    hs_rs_with_tip = @view hs_rs[:, end:end, :]
    Psis_rs_with_tip = @view Psis_rs[:, end:end, :]
    bls_rs_with_tip = @view bls_rs[:, end:end, :]
    trans_with_tip = @view trans[:, end:end]

    # Construct and transform the source elements.
    ses_no_tip = CombinedNoTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}.(Ref(rotor), sections_rs_no_tip, ops_rs_no_tip, outputs_rs_no_tip, θs_rs, dradii_rs_no_tip, hs_rs_no_tip, Psis_rs_no_tip, src_times, Ref(dt), bls_rs_no_tip, positive_x_rotation) .|> trans_no_tip

    ses_with_tip = CombinedWithTipBroadbandSourceElement{TDirect,TUInduction,TMachCorrection,TDoppler}.(Ref(rotor), sections_rs_with_tip, ops_rs_with_tip, outputs_rs_with_tip, θs_rs, dradii_rs_with_tip, hs_rs_with_tip, Psis_rs_with_tip, src_times, Ref(dt), bls_rs_with_tip, Ref(blade_tip), positive_x_rotation) .|> trans_with_tip

    return ses_no_tip, ses_with_tip
end

"""
    get_ccblade_dradii(rotor::CCBlade.Rotor, sections::Vector{CCBlade.Section})

Construct and return a Vector of the lengths of each CCBlade section.
"""
function get_ccblade_dradii(rotor, sections)
    radii = mapview(:r, sections)
    dradii = get_dradii(radii, rotor.Rhub, rotor.Rtip)
    return dradii
end
