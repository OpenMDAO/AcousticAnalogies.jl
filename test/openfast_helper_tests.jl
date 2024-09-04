module OpenFASTHelperTests

using SafeTestsets: @safetestset
using Test: @testset


@testset "OpenFAST reader tests" begin

    @testset "reading" begin
        @safetestset "default" begin
            include("openfast_reading_default.jl")
        end

        @safetestset "different time column name" begin
            include("openfast_reading_different_time_column_name.jl")
        end

        @safetestset "no units header" begin
            include("openfast_reading_no_units_header.jl")
        end
    end

    @testset "radial interpolation" begin
        @safetestset "linear" begin
            include("openfast_radial_interpolation_linear.jl")
        end

        @safetestset "akima" begin
            include("openfast_radial_interpolation_akima.jl")
        end
    end

    @testset "time derivatives" begin
        @safetestset "constant time step" begin
            include("openfast_time_derivatives_constant_time_step.jl")
        end

        @safetestset "non-constant time step" begin
            include("openfast_time_derivatives_nonconstant_time_step.jl")
        end
        
        @safetestset "no time diff" begin
            include("openfast_time_derivatives_nonconstant_time_step_no_diff.jl")
        end
    end

end

end
