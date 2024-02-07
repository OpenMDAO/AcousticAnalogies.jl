abstract type AbstractCompactSourceElement end

"""
    velocity(se::AbstractCompactSourceElement)

Return the current velocity of `se`.
"""
@inline velocity(se::AbstractCompactSourceElement) = se.y1dot

