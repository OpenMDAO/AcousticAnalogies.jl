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
