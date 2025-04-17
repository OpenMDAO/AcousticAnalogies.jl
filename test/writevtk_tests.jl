module WriteVTKTests

using AcousticAnalogies
using Format: format, FormatExpr
using JLD2: JLD2
using SHA: sha1
using StaticArrays: @SVector
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
            fname = format(FormatExpr("{}{:08d}.vtp"), name, i)
            sha_str = bytes2hex(open(sha1, fname))
            sha_str_check = bytes2hex(open(sha1, joinpath("writevtk", fname)))
            @test sha_str == sha_str_check
        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        if !Sys.iswindows()
            @test sha_str == sha_str_check
        end

    end

    @testset "Compact F1A source elements, with observers" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            ses = file["ses"]
        end

        obs1 = AcousticAnalogies.ConstVelocityAcousticObserver(0.0, @SVector([0, 2.0, 0]), @SVector([5.0, 0.0, 0.0]))
        obs2 = AcousticAnalogies.StationaryAcousticObserver(@SVector [0, 2.5, 0])
        obs = [obs1, obs2]

        name = "cf1a_with_observers"
        pvd = AcousticAnalogies.to_paraview_collection(name, (ses,); observers=obs, append=false, ascii=true)

        for i in 1:size(ses, 1)
            fname = format(FormatExpr("{}-block1-{:08d}.vtp"), name, i)
            sha_str = bytes2hex(open(sha1, fname))
            sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
            if !Sys.iswindows()
                @test sha_str == sha_str_check
            end

            # The source element files for this test case with observers should be the same as the case without the observers.
            name2 = "cf1a"
            fname2 = format(FormatExpr("{}{:08d}.vtp"), name2, i)
            sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
            if !Sys.iswindows()
                @test sha_str == sha_str_check
            end

            for j in 1:length(obs)
                # This just isn't stable across versions of Meshes.jl.
                # I suppose it has something to do with Delaunay triangulation.
                fname = format(FormatExpr("{}-observer$(j)-{:08d}.vtu"), name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                # @test sha_str == sha_str_check
            end

        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        if !Sys.iswindows()
            @test sha_str == sha_str_check
        end

    end

    @testset "Compact F1A source elements, multiblock, with observers" begin
        fname = joinpath(@__DIR__, "writevtk", "cf1a.jld2")

        ses = nothing
        JLD2.jldopen(fname, "r") do file
            ses = file["ses"]
        end

        # Split the array into "blocks."
        ses_mb = tuple([ses[:, :, b] for b in 1:size(ses, 3)]...)

        obs1 = AcousticAnalogies.ConstVelocityAcousticObserver(0.0, @SVector([0, 2.0, 0]), @SVector([5.0, 0.0, 0.0]))
        obs2 = AcousticAnalogies.StationaryAcousticObserver(@SVector [0, 2.5, 0])
        obs = [obs1, obs2]

        name = "cf1a_mb_with_observers"
        pvd = AcousticAnalogies.to_paraview_collection(name, ses_mb; observers=obs, append=false, ascii=true)

        for i in 1:size(ses, 1)
            for b in 1:length(ses_mb)
                fname = format(FormatExpr("{}-block$(b)-{:08d}.vtp"), name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                if !Sys.iswindows()
                    @test sha_str == sha_str_check
                end
            end

            for j in 1:length(obs)
                # This just isn't stable across versions of Meshes.jl.
                # I suppose it has something to do with Delaunay triangulation.
                fname = format(FormatExpr("{}-observer$(j)-{:08d}.vtu"), name, i)
                sha_str = bytes2hex(open(sha1, fname))
                sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
                # @test sha_str == sha_str_check
            end

        end

        fname = "$(name).pvd"
        sha_str = bytes2hex(open(sha1, fname))
        sha_str_check = bytes2hex(open(sha1, joinpath(@__DIR__, "writevtk", fname)))
        if !Sys.iswindows()
            @test sha_str == sha_str_check
        end

    end

end

end # module
