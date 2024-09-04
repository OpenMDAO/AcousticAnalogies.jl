module AcousticAnalogiesTests

all_tests = isempty(ARGS) || ("all" in ARGS)

if "adv_time" in ARGS || all_tests
    include("adv_time_tests.jl")
end
if "combine" in ARGS || all_tests
    include("combine_tests.jl")
end
if "f1a" in ARGS || all_tests
    include("f1a_tests.jl")
end
if "f1a_constructor" in ARGS || all_tests
    include("compact_f1a_constructor_tests.jl")
end
if "anopp2" in ARGS || all_tests
    include("anopp2_comparison.jl")
end
if "forwarddiff" in ARGS || all_tests
    include("forwarddiff_test.jl")
end
if "doppler" in ARGS || all_tests
    include("doppler_tests.jl")
end
if "boundary_layers" in ARGS || all_tests
    include("boundary_layer_tests.jl")
end
if "bpm_shape_functions" in ARGS || all_tests
    include("bpm_shape_function_tests.jl")
end
if "broadband_source_elements" in ARGS || all_tests
    include("broadband_source_element_tests.jl")
end
if "writevtk" in ARGS || all_tests
    include("writevtk_tests.jl")
end
if "openfast" in ARGS || all_tests
    include("openfast_helper_tests.jl")
end
if "bpm_itr" in ARGS || all_tests
    include("bpm_itr_tests.jl")
end

end # module
