module WriteVTKTests

using AcousticAnalogies
using Formatting: format
using JLD2: JLD2
using SHA: sha1
using Test

@testset "WriteVTK tests" begin

    @testset "Compact F1A source elements" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            ses = file["ses"]
        end

        name = "cf1a"
        pvd = AcousticAnalogies.to_paraview_collection(name, ses)

        for i in 1:size(ses, 1)
            fname = format("{}{:08d}.vtp", name, i)
            sha_str = bytes2hex(open(sha1, fname))
            sha_str_check = bytes2hex(open(sha1, joinpath("writevtk", fname)))
            @test sha_str == sha_str_check
        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath("writevtk", fname)))
        @test sha_str == sha_str_check

    end

end

end # module
