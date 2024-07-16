module AcousticAnalogies

using AcousticMetrics: AcousticMetrics
using BPM: BPM
using CCBlade: CCBlade
using FillArrays: Fill
using FLOWMath: akima, linear, ksmax, norm_cs_safe, dot_cs_safe, atan_cs_safe, abs_cs_safe
using FlexiMaps: mapview
using Formatting: format
using JuliennedArrays: JuliennedArrays
using KinematicCoordinateTransformations: KinematicTransformation, SteadyRotXTransformation, ConstantVelocityTransformation, compose
using LinearAlgebra: cross, norm, mul!
using Meshes: Meshes
using StaticArrays: @SVector
using WriteVTK: WriteVTK

include("utils.jl")
export get_dradii

include("abstract_source_elements.jl")
export AbstractCompactSourceElement

include("observers.jl")
export AbstractAcousticObserver, StationaryAcousticObserver, ConstVelocityAcousticObserver

include("advance_time.jl")
export adv_time

include("boundary_layers.jl")
export AbstractBoundaryLayer, TrippedN0012BoundaryLayer, UntrippedN0012BoundaryLayer

include("f1a.jl")
export CompactF1ASourceElement
export F1AOutput, F1APressureTimeHistory
export noise
export common_obs_time
export combine!, combine

include("abstract_broadband.jl")

include("tbl_te.jl")
export TBLTESourceElement

include("lbl_vs.jl")
export LBLVSSourceElement

include("tip_vortex.jl")
export AbstractTipAlphaCorrection, NoTipAlphaCorrection, BPMTipAlphaCorrection, BMTipAlphaCorrection, SmoothBMTipAlphaCorrection
export AbstractBladeTip, RoundedTip, FlatTip
export TipVortexSourceElement

include("teb_vs.jl")
export TEBVSSourceElement

include("combined_broadband.jl")
export CombinedNoTipBroadbandSourceElement, CombinedWithTipBroadbandSourceElement
export pbs_suction, pbs_pressure, pbs_alpha, pbs_teb, pbs_tip

include("ccblade_helpers.jl")
export f1a_source_elements_ccblade, tblte_source_elements_ccblade, lblvs_source_elements_ccblade, tebvs_source_elements_ccblade, tip_vortex_source_elements_ccblade, combined_broadband_source_elements_ccblade

include("bpm_test_utils.jl")

include("writevtk.jl")

end # module
