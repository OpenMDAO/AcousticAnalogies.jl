module BPMITRTests

using SafeTestsets: @safetestset
using Test: @testset

@testset "BPM.jl comparisons" begin
    @safetestset "figure 22b" begin
        include("itr_figure22b_bpmjl.jl")
    end

    @safetestset "figure 23c" begin
        include("itr_figure23c_bpmjl.jl")
    end

    @safetestset "figure 24b" begin
        include("itr_figure24b_bpmjl.jl")
    end

end

@testset "PAS/ROTONET/BARC comparisons" begin
    @safetestset "figure 22b" begin
        include("itr_figure22b_barc.jl")
    end

    @safetestset "figure 23c" begin
        include("itr_figure23c_barc.jl")
    end

    @safetestset "figure 24b" begin
        include("itr_figure24b_barc.jl")
    end
end

end # module
