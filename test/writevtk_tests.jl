module WriteVTKTests

using AcousticAnalogies
using Formatting: format
using JLD2: JLD2
using SHA: sha1
using StaticArrays: @SVector
using Test

function to_cf1a(se)
    return CompactF1ASourceElement(se.ρ0, se.c0, se.Δr, se.Λ, se.y0dot, se.y1dot, se.y2dot, se.y3dot, se.f0dot, se.f1dot, se.τ, se.u)
end

@testset "WriteVTK tests" begin

    @testset "Compact F1A source elements" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            # Renaming CompactSourceElement to CompactF1ASourceElement breaks reconstructing the source elements from the jld2file.
            ses = to_cf1a.(file["ses"])
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
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        @test sha_str == sha_str_check

    end

    @testset "Compact F1A source elements, with observers" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            # Renaming CompactSourceElement to CompactF1ASourceElement breaks reconstructing the source elements from the jld2file.
            ses = to_cf1a.(file["ses"])
        end

        obs1 = AcousticAnalogies.ConstVelocityAcousticObserver(0.0, @SVector([0, 2.0, 0]), @SVector([5.0, 0.0, 0.0]))
        obs2 = AcousticAnalogies.StationaryAcousticObserver(@SVector [0, 2.5, 0])
        obs = [obs1, obs2]

        name = "cf1a_with_observers"
        pvd = AcousticAnalogies.to_paraview_collection(name, (ses,); observers=obs)

        for i in 1:size(ses, 1)
            fname = format("{}-block1-{:08d}.vtp", name, i)
            sha_str = bytes2hex(open(sha1, fname))
            sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
            @test sha_str == sha_str_check

            # The source element files for this test case with observers should be the same as the case without the observers.
            name2 = "cf1a"
            fname2 = format("{}{:08d}.vtp", name2, i)
            sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
            @test sha_str == sha_str_check

            for j in 1:length(obs)
                fname = format("{}-observer$(j)-{:08d}.vtu", name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                # @test sha_str == sha_str_check
            end

        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        @test sha_str == sha_str_check

    end

    @testset "Compact F1A source elements, multiblock, with observers" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            # Renaming CompactSourceElement to CompactF1ASourceElement breaks reconstructing the source elements from the jld2file.
            ses = to_cf1a.(file["ses"])
        end

        # Split the array into "blocks."
        ses_mb = tuple([ses[:, :, b] for b in 1:size(ses, 3)]...)

        obs1 = AcousticAnalogies.ConstVelocityAcousticObserver(0.0, @SVector([0, 2.0, 0]), @SVector([5.0, 0.0, 0.0]))
        obs2 = AcousticAnalogies.StationaryAcousticObserver(@SVector [0, 2.5, 0])
        obs = [obs1, obs2]

        name = "cf1a_mb_with_observers"
        pvd = AcousticAnalogies.to_paraview_collection(name, ses_mb; observers=obs)

        for i in 1:size(ses, 1)
            for b in 1:length(ses_mb)
                fname = format("{}-block$(b)-{:08d}.vtp", name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                @test sha_str == sha_str_check
            end

            for j in 1:length(obs)
                fname = format("{}-observer$(j)-{:08d}.vtu", name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                # @test sha_str == sha_str_check

                # The observers for this case should be identical to the observers from the single-block case.
                fname2 = format("cf1a_with_observers-observer$(j)-{:08d}.vtu", i)
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname2)))
                # @test sha_str == sha_str_check
            end

        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        @test sha_str == sha_str_check

    end

end

end # module
