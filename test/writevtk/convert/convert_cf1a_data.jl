module ConvertCF1AData

using AcousticAnalogies
using JLD2: JLD2

function to_cf1a(se)
    return CompactF1ASourceElement(se.ρ0, se.c0, se.Δr, se.Λ, se.y0dot, se.y1dot, se.y2dot, se.y3dot, se.f0dot, se.f1dot, se.τ, se.u)
end

function doit()
    ses = nothing
    JLD2.jldopen(joinpath(@__DIR__, "cf1a-old.jld2"), "r") do file
        # Renaming CompactSourceElement to CompactF1ASourceElement breaks reconstructing the source elements from the jld2file.
        ses = to_cf1a.(file["ses"])
    end

    JLD2.jldopen(joinpath(@__DIR__, "cf1a.jld2"), "w") do file
        file["ses"] = ses
    end

end

end # module
