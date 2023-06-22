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

abstract type AbstractCompactSourceElement end

include("boundary_layers.jl")
export AbstractBoundaryLayer, TrippedN0012BoundaryLayer, UntrippedN0012BoundaryLayer

include("core.jl")
export CompactSourceElement
export AcousticObserver, StationaryAcousticObserver, ConstVelocityAcousticObserver
export F1AOutput, F1APressureTimeHistory
export adv_time
export f1a
export common_obs_time
export combine!, combine

include("tbl_te.jl")
include("lbl_vs.jl")
include("tip_vortex.jl")
include("teb_vs.jl")

include("ccblade_helpers.jl")
export source_elements_ccblade, tblte_source_elements_ccblade

include("writevtk.jl")

end # module
