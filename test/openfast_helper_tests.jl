module OpenFASTHelperTests

using AcousticAnalogies
using Statistics: mean
using Test

@testset "OpenFAST reader test" begin

    @testset "default" begin
        fname = joinpath(@__DIR__, "gen_test_data", "openfast_data", "IEA-3.4-130-RWT-small.out")

        data = read_openfast_file(fname)
        num_times = length(data.time)
        num_blades = size(data.pitch, 2)
        num_radial = size(data.axial_loading, 2)
        @test size(data.time) == (num_times,)
        @test size(data.v) == (num_times,)
        @test size(data.azimuth) == (num_times,)
        @test size(data.omega) == (num_times,)
        @test size(data.pitch) == (num_times, num_blades)
        @test size(data.axial_loading) == (num_times, num_radial, num_blades)
        @test size(data.circum_loading) == (num_times, num_radial, num_blades)

        @test all(data.time .≈ (60.0:0.01:60.1))
        @test all(data.v .≈ 7)
        @test all(data.pitch .≈ 0)
        # Just some coarse tests.
        @test mean(data.omega) ≈ 8.126818181818182
        @test mean(data.axial_loading) ≈ 1632.1274949494953
        @test mean(data.circum_loading) ≈ -217.45748949494939
    end

    @testset "different time column name" begin
        fname = joinpath(@__DIR__, "gen_test_data", "openfast_data", "IEA-3.4-130-RWT-small-FooTime.out")

        data = read_openfast_file(fname; header_keyword="FooTime")
        num_times = length(data.time)
        num_blades = size(data.pitch, 2)
        num_radial = size(data.axial_loading, 2)
        @test size(data.time) == (num_times,)
        @test size(data.v) == (num_times,)
        @test size(data.azimuth) == (num_times,)
        @test size(data.omega) == (num_times,)
        @test size(data.pitch) == (num_times, num_blades)
        @test size(data.axial_loading) == (num_times, num_radial, num_blades)
        @test size(data.circum_loading) == (num_times, num_radial, num_blades)

        @test all(data.time .≈ (60.0:0.01:60.1))
        @test all(data.v .≈ 7)
        @test all(data.pitch .≈ 0)
        # Just some coarse tests.
        @test mean(data.omega) ≈ 8.126818181818182
        @test mean(data.axial_loading) ≈ 1632.1274949494953
        @test mean(data.circum_loading) ≈ -217.45748949494939
    end

    @testset "nu units header" begin
        fname = joinpath(@__DIR__, "gen_test_data", "openfast_data", "IEA-3.4-130-RWT-small-no_units.out")

        data = read_openfast_file(fname; has_units_header=false)
        num_times = length(data.time)
        num_blades = size(data.pitch, 2)
        num_radial = size(data.axial_loading, 2)
        @test size(data.time) == (num_times,)
        @test size(data.v) == (num_times,)
        @test size(data.azimuth) == (num_times,)
        @test size(data.omega) == (num_times,)
        @test size(data.pitch) == (num_times, num_blades)
        @test size(data.axial_loading) == (num_times, num_radial, num_blades)
        @test size(data.circum_loading) == (num_times, num_radial, num_blades)

        @test all(data.time .≈ (60.0:0.01:60.1))
        @test all(data.v .≈ 7)
        @test all(data.pitch .≈ 0)
        # Just some coarse tests.
        @test mean(data.omega) ≈ 8.126818181818182
        @test mean(data.axial_loading) ≈ 1632.1274949494953
        @test mean(data.circum_loading) ≈ -217.45748949494939
    end

end

end
