abstract type AbstractCompactSourceElement end

@inline velocity(se::AbstractCompactSourceElement) = se.y1dot

