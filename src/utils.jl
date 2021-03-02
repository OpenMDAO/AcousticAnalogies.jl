function norm_cs_safe(x::AbstractVector{Complex{T}}) where {T}
    return sqrt(sum(x.*x))
end

function norm_cs_safe(x)
    return norm(x)
end

function dot_cs_safe(a, b)
    return sum(a.*b)
end
