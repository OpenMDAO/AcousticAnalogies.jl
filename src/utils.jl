# Smooth maximum: https://en.wikipedia.org/wiki/LogSumExp,
# https://www.johndcook.com/blog/2010/01/13/soft-maximum/
function smooth_max(x::AbstractArray{T}; k=one(T)) where {T}
    return log(sum(exp.(k*x)))/k
end

function smooth_max_time(apth, k, t0)
    return log(sum(exp.(k.*(getproperty.(apth, :t).-t0))))/k
end

function norm_cs_safe(x::AbstractVector{Complex{T}}) where {T}
    # return norm(real(x)) + norm(imag(x))im
    return sqrt(sum(x.*x))
end

function norm_cs_safe(x)
    return norm(x)
end

function dot_cs_safe(a, b)
    return sum(a.*b)
end
