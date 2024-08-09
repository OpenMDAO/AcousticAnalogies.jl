abstract type AbstractCompactSourceElement end

"""
    velocity(se::AbstractCompactSourceElement)

Return the current velocity of `se`.
"""
@inline velocity(se::AbstractCompactSourceElement) = se.y1dot

"""
    source_time(se::AbstractCompactSourceElement)

Return the source time of `se`.
"""
@inline source_time(se::AbstractCompactSourceElement) = se.τ

"""
    orientation(se::AbstractCompactSourceElement)

Return a length-3 unit vector indicating the spanwise orientation of `se`.
"""
@inline orientation(se::AbstractCompactSourceElement) = se.span_uvec

"""
    position(se::AbstractCompactSourceElement)

Return a length-3 vector indicating the position of `se`.
"""
@inline position(se::AbstractCompactSourceElement) = se.y0dot

"""
    speed_of_sound(se::AbstractCompactSourceElement)

Return the ambient speed of sound associated with `se`.
"""
@inline speed_of_sound(se) = se.c0
