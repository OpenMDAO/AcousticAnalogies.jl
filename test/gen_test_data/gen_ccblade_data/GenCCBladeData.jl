module GenCCBladeData

using CCBlade: CCBlade
using FileIO: save
using JLD2: JLD2
using Printf: @sprintf

include("constants.jl")
include("xrotor_airfoil.jl")

# Most of this is taken from <developer.nasa.gov:dingraha/CROTORCCBladeComparisons.jl.git>

#' This is the function that will iterate over the RPM values and do the computation.
function omega_sweep()
    #' I ended up extracting and CROTOR airfoil code and using it with CCBlade to make
    #' sure we're actually doing an "apples-to-apples" comparsion. This is a Julia
    #' `struct` that holds all the parameters needed for the CROTOR airfoil routine.
    #' The `af_xrotor` routine takes the `XROTORAirfoilConfig` `struct` and
    #' an angle of attack, Reynolds number, and Mach number and returns a lift and
    #' drag coefficient.
    xrotor_config = XROTORAirfoilConfig(
        A0=0.0, DCLDA=6.2800, CLMAX=1.5, CLMIN=-0.5, DCLDA_STALL=0.1,
        DCL_STALL=0.1, MCRIT=0.8, CDMIN=0.13e-1, CLDMIN=0.5, DCDCL2=0.4e-2, REREF=0.2e6, REXP=-0.4)

    airfoil_interp(a, r, m) = af_xrotor(a, r, m, xrotor_config)
    theta_rad = CCBladeTestCaseConstants.theta .* pi/180.0

    rotor = CCBlade.Rotor(CCBladeTestCaseConstants.Rhub, CCBladeTestCaseConstants.Rtip, CCBladeTestCaseConstants.num_blades)
    sections = CCBlade.Section.(CCBladeTestCaseConstants.radii, CCBladeTestCaseConstants.chord, theta_rad, airfoil_interp)

    rpm = 200.0:200.0:2200.0  # rev/min
    for i in eachindex(rpm)
        omega = rpm[i]*(2*pi/60.0)

        ops = CCBlade.OperatingPoint.(CCBladeTestCaseConstants.v, omega.*CCBladeTestCaseConstants.radii, CCBladeTestCaseConstants.rho, CCBladeTestCaseConstants.pitch, CCBladeTestCaseConstants.mu, CCBladeTestCaseConstants.c0)
        outs = CCBlade.solve.(Ref(rotor), sections, ops)

        outs_d = Dict{String,Any}(string(f)=>getproperty.(outs, f) for f in fieldnames(eltype(outs)))
        outs_d["omega"] = omega
        save("ccblade_omega$(@sprintf "%02d" i)-outputs.jld2", outs_d)

        # data = DelimitedFiles.readdlm("ccblade_omega$(@sprintf "%02d" i).csv", ',')
        # @show maximum(abs.(getproperty.(outs, :Np) .- data[:, 1]))
        # @show maximum(abs.(getproperty.(outs, :Tp) .- data[:, 2]))
    end 

    return nothing
end

if ! isinteractive()
    omega_sweep()
end

end # module
