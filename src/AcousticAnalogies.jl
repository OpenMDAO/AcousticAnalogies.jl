module AcousticAnalogies

using AcousticMetrics
using CCBlade
using ConcreteStructs: @concrete
using FillArrays: Fill
using FLOWMath: akima, linear, ksmax, norm_cs_safe, dot_cs_safe, atan_cs_safe, abs_cs_safe
using Formatting: format
using KinematicCoordinateTransformations
using LinearAlgebra: cross, norm, mul!
using SingleFieldStructArrays
using StaticArrays
using WriteVTK

include("utils.jl")
export get_dradii

include("abstract_source_elements.jl")

include("boundary_layers.jl")
export AbstractBoundaryLayer, TrippedN0012BoundaryLayer, UntrippedN0012BoundaryLayer

include("core.jl")
export CompactSourceElement
export AbstractAcousticObserver, StationaryAcousticObserver, ConstVelocityAcousticObserver
export F1AOutput, F1APressureTimeHistory
export adv_time
export f1a
export common_obs_time
export combine!, combine

include("tbl_te.jl")
export TBLTESourceElement
include("lbl_vs.jl")
export LBLVSSourceElement
include("tip_vortex.jl")
export TipVortexSourceElement
include("teb_vs.jl")
export TEBVSSourceElement
include("combined_broadband.jl")
export CombinedNoTipBroadbandSourceElement, CombinedWithTipBroadbandSourceElement
export pbs_suction, pbs_pressure, pbs_alpha, pbs_teb, pbs_tip

include("ccblade_helpers.jl")
export source_elements_ccblade, tblte_source_elements_ccblade, lblvs_source_elements_ccblade, tebvs_source_elements_ccblade, tip_vortex_source_elements_ccblade, combined_broadband_source_elements_ccblade

include("bpm_test_utils.jl")

include("writevtk.jl")

end # module
