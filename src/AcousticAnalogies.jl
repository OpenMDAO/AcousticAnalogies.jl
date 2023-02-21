module AcousticAnalogies

using AcousticMetrics
using CCBlade
using ConcreteStructs: @concrete
using FLOWMath: akima, linear, ksmax
using Formatting: format
using KinematicCoordinateTransformations
using LinearAlgebra: cross, norm, mul!
using SingleFieldStructArrays
using StaticArrays
using WriteVTK

include("utils.jl")
export get_dradii

include("boundary_layers.jl")
include("shape_functions.jl")

include("core.jl")
export CompactSourceElement
export AcousticObserver, StationaryAcousticObserver, ConstVelocityAcousticObserver
export F1AOutput, F1APressureTimeHistory
export adv_time
export f1a
export common_obs_time
export combine!, combine

include("ccblade_helpers.jl")
export source_elements_ccblade

include("writevtk.jl")

end # module
