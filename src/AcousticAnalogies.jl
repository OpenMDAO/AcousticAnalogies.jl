module AcousticAnalogies

using AcousticMetrics
using CCBlade
using ConcreteStructs: @concrete
using FLOWMath: akima, linear, ksmax, dot_cs_safe, norm_cs_safe
using Formatting: format
using JuliennedArrays: JuliennedArrays
using KinematicCoordinateTransformations
using LinearAlgebra: cross, norm, mul!
using FlexiMaps: mapview
using StaticArrays
using WriteVTK

include("utils.jl")
export get_dradii

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
