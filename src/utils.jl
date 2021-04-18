function norm_cs_safe(x::AbstractVector{Complex{T}}) where {T}
    return sqrt(sum(x.*x))
end

function norm_cs_safe(x)
    return norm(x)
end

function dot_cs_safe(a, b)
    return sum(a.*b)
end

"""
    get_dradii(radii, Rhub, Rtip)

Compute the spacing between blade elements given the radial locations of the
element midpoints in `radii` and the hub and tip radius in `Rhub` and `Rtip`,
respectively.

Assume the interfaces between elements are midway between adjacent element's midpoints.
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
