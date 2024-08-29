abstract type AbstractTimeDerivMethod end
struct NoTimeDerivMethod <: AbstractTimeDerivMethod end
struct SecondOrderFiniteDiff <: AbstractTimeDerivMethod end

abstract type AbstractRadialInterpMethod end
struct FLOWLinearInterp <: AbstractRadialInterpMethod end
struct FLOWAkimaInterp <: AbstractRadialInterpMethod end

"""
Struct for holding data from an OpenFAST (AeroDyn?) output file.

# Fields:
* `time`: vector of simulation times with size `(num_times,)`
* `dtime_dtau`: vector of derivative of simulation times with respect to non-dimensional/computational time with size `(num_times,)`
* `v`: vector of freestream velocity time history with size `(num_times,)`
* `azimuth`: vector of azimuth angle time history with size `(num_times,)`
* `omega`: vector of rotation rate time history with size `(num_times,)`
* `pitch`: array of pitch angle time history with size `(num_times, num_blades)`
* `radii`: vector of blade radial locations with size `(num_radial,)`
* `radii_mid`: vector of cell-centered/midpoint blade radial locations with size `(num_radial-1,)`
* `cs_area`: vector of cross-sectional areas with size `(num_radial,)`.
* `cs_area_mid`: vector of cell-centered/midpoint cross-sectional areas with size `(num_radial-1,)`.
* `axial_loading`: array of axial loading time history with size `(num_times, num_radial, num_blades)`
* `axial_loading_mid`: array of axial loading time history at cell-centered blade radial locations with size `(num_times, num_radial-1, num_blades)`
* `axial_loading_mid_dot`: array of axial loading temporal derivative time history at cell-centered blade radial locations with size `(num_times, num_radial-1, num_blades)`
* `circum_loading`: array of circumferential loading time history with size `(num_times, num_radial, num_blades)`
* `circum_loading_mid`: array of circum loading time history at cell-centered blade radial locations with size `(num_times, num_radial-1, num_blades)`
* `circum_loading_mid_dot`: array of circum loading temporal derivative time history at cell-centered blade radial locations with size `(num_times, num_radial-1, num_blades)`
"""
struct OpenFASTData{TRadialInterpMethod,TTimeDerivMethod,
                    TTime,TdTimedTau,TV,TAzimuth,TOmega,TPitch,
                    TRadii,TRadiiMid,TDRadii,
                    TCSArea,TCSAreaMid,
                    TAxialLoading,TAxialLoadingMid,TAxialLoadingMidDot,
                    TCircumLoading,TCircumLoadingMid,TCircumLoadingMidDot}
    time::TTime
    dtime_dtau::TdTimedTau
    v::TV
    azimuth::TAzimuth
    omega::TOmega
    pitch::TPitch
    radii::TRadii
    radii_mid::TRadiiMid
    dradii::TDRadii
    cs_area::TCSArea
    cs_area_mid::TCSAreaMid
    axial_loading::TAxialLoading
    axial_loading_mid::TAxialLoadingMid
    axial_loading_mid_dot::TAxialLoadingMidDot
    circum_loading::TCircumLoading
    circum_loading_mid::TCircumLoadingMid
    circum_loading_mid_dot::TCircumLoadingMidDot
    num_blades::Int

    function OpenFASTData{
        TRadialInterpMethod,TTimeDerivMethod,
        TTime,TdTimedTau,TV,TAzimuth,TOmega,TPitch,
        TRadii,TRadiiMid,TDRadii,
        TCSArea,TCSAreaMid,
        TAxialLoading,TAxialLoadingMid,TAxialLoadingMidDot,
        TCircumLoading,TCircumLoadingMid,TCircumLoadingMidDot}(
            time, dtime_dtau, v, azimuth, omega, pitch,
            radii, radii_mid, dradii,
            cs_area, cs_area_mid,
            axial_loading, axial_loading_mid, axial_loading_mid_dot,
            circum_loading, circum_loading_mid, circum_loading_mid_dot) where {
                TRadialInterpMethod,TTimeDerivMethod,
                TTime,TdTimedTau,TV,TAzimuth,TOmega,TPitch,
                TRadii,TRadiiMid,TDRadii,
                TCSArea,TCSAreaMid,
                TAxialLoading,TAxialLoadingMid,TAxialLoadingMidDot,
                TCircumLoading,TCircumLoadingMid,TCircumLoadingMidDot}

        num_times = length(time)
        num_radial = length(radii)
        num_radial_mid = num_radial - 1

        # Figure out what num_blades is.
        if pitch !== nothing
            num_blades = size(pitch, 2)
        elseif axial_loading !== nothing
            num_blades = size(axial_loading, 3)
        elseif axial_loading_mid !== nothing
            num_blades = size(axial_loading_mid, 3)
        elseif axial_loading_mid_dot !== nothing
            num_blades = size(axial_loading_mid_dot, 3)
        elseif circum_loading !== nothing
            num_blades = size(circum_loading, 3)
        elseif circum_loading_mid !== nothing
            num_blades = size(circum_loading_mid, 3)
        elseif circum_loading_mid_dot !== nothing
            num_blades = size(circum_loading_mid_dot, 3)
        else
            num_blades = 0
        end

        if dtime_dtau !== nothing
            size(dtime_dtau) == (num_times,) || throw(ArgumentError("size(dtime_dtau) = $(size(dtime_dtau)) does not match size(time) = $(size(time))"))
        end

        if v !== nothing
            size(v) == (num_times,) || throw(ArgumentError("size(v) = $(size(v)) does not match size(time) = $(size(time))"))
        end

        if azimuth !== nothing
            size(azimuth) == (num_times,) || throw(ArgumentError("size(azimuth) = $(size(azimuth)) does not match size(time) = $(size(time))"))
        end

        if omega !== nothing
            size(omega) == (num_times,) || throw(ArgumentError("size(omega) = $(size(omega)) does not match size(time) = $(size(time))"))
        end

        if pitch !== nothing
            size(pitch) == (num_times, num_blades) || throw(ArgumentError("size(pitch) = $(size(pitch)) not consistent with size(time) = $(size(time)) and/or size(axial_loading) = $(size(axial_loading))"))
        end

        if radii_mid !== nothing
            size(radii_mid) == (num_radial_mid,) || throw(ArgumentError("size(radii_mid) = $(size(radii_mid)) not consistent with length(radii)-1 = $(num_radial_mid)"))
        end

        if dradii !== nothing
            size(dradii) == (num_radial_mid,) || throw(ArgumentError("size(dradii) = $(size(dradii)) not consistent with length(radii)-1 = $(num_radial_mid)"))
        end

        if cs_area !== nothing
            size(cs_area) == (num_radial,) || throw(ArgumentError("size(cs_area) = $(size(cs_area)) not consistent with length(radii) = $(num_radial)"))
        end

        if cs_area_mid !== nothing
            size(cs_area_mid) == (num_radial_mid,) || throw(ArgumentError("size(cs_area_mid) = $(size(cs_area_mid)) not consistent with length(radii)-1 = $(num_radial_mid)"))
        end

        if axial_loading !== nothing
            size(axial_loading) == (num_times, num_radial, num_blades) || throw(ArgumentError("size(axial_loading) = $(size(axial_loading)) not consistent with size(time) = $(size(time)), length(radii) = $(num_radial), blade count = $(num_blades)"))
        end

        if axial_loading_mid !== nothing
            size(axial_loading_mid) == (num_times, num_radial_mid, num_blades) || throw(ArgumentError("size(axial_loading_mid) = $(size(axial_loading_mid)) not consistent with size(time) = $(size(time)), length(radii)-1 = $(num_radial_mid), blade count = $(num_blades)"))
        end

        if axial_loading_mid_dot !== nothing
            size(axial_loading_mid_dot) == (num_times, num_radial_mid, num_blades) || throw(ArgumentError("size(axial_loading_mid_dot) = $(size(axial_loading_mid_dot)) not consistent with size(time) = $(size(time)), length(radii)-1 = $(num_radial_mid), blade count = $(num_blades)"))
        end

        if circum_loading !== nothing
            size(circum_loading) == (num_times, num_radial, num_blades) || throw(ArgumentError("size(circum_loading) = $(size(circum_loading)) not consistent with size(time) = $(size(time)), length(radii) = $(num_radial), blade count = $(num_blades)"))
        end

        if circum_loading_mid !== nothing
            size(circum_loading_mid) == (num_times, num_radial_mid, num_blades) || throw(ArgumentError("size(circum_loading_mid) = $(size(circum_loading_mid)) not consistent with size(time) = $(size(time)), length(radii)-1 = $(num_radial_mid), blade count = $(num_blades)"))
        end

        if circum_loading_mid_dot !== nothing
            size(circum_loading_mid_dot) == (num_times, num_radial_mid, num_blades) || throw(ArgumentError("size(circum_loading_mid_dot) = $(size(circum_loading_mid_dot)) not consistent with size(time) = $(size(time)), length(radii)-1 = $(num_radial_mid), blade count = $(num_blades)"))
        end

        return new{TRadialInterpMethod,TTimeDerivMethod,
                   typeof(time),typeof(dtime_dtau),typeof(v),typeof(azimuth),typeof(omega),typeof(pitch),
                   typeof(radii),typeof(radii_mid),typeof(dradii),
                   typeof(cs_area),typeof(cs_area_mid),
                   typeof(axial_loading),typeof(axial_loading_mid),typeof(axial_loading_mid_dot),
                   typeof(circum_loading),typeof(circum_loading_mid),typeof(circum_loading_mid_dot)}(
                       time, dtime_dtau, v, azimuth, omega, pitch,
                       radii, radii_mid, dradii,
                       cs_area, cs_area_mid,
                       axial_loading, axial_loading_mid, axial_loading_mid_dot,
                       circum_loading, circum_loading_mid, circum_loading_mid_dot,
                      num_blades)
    end
end

function OpenFASTData{TRadialInterpMethod,TTimeDerivMethod}(time, v, azimuth, omega, pitch, radii, cs_area, axial_loading, circum_loading) where {TRadialInterpMethod<:AbstractRadialInterpMethod,TTimeDerivMethod<:AbstractTimeDerivMethod}
    dtime_dtau = similar(time)

    # Find the blade element midpoint locations.
    radii_mid = 0.5 .* (@view(radii[begin:end-1]) .+ @view(radii[begin+1:end]))

    # Find the radial spacing.
    dradii = @view(radii[begin+1:end]) .- @view(radii[begin:end-1])

    if cs_area !== nothing
        cs_area_mid = similar(cs_area, length(cs_area)-1)
    else
        cs_area_mid = nothing
    end

    if axial_loading !== nothing
        axial_loading_mid = similar(axial_loading, size(axial_loading, 1), size(axial_loading, 2)-1, size(axial_loading, 3))
        axial_loading_mid_dot = zero(axial_loading_mid)
    else
        axial_loading_mid = axial_loading_mid_dot = nothing
    end

    if circum_loading !== nothing
        circum_loading_mid = similar(circum_loading, size(circum_loading, 1), size(circum_loading, 2)-1, size(circum_loading, 3))
        circum_loading_mid_dot = zero(circum_loading_mid)
    else
        circum_loading_mid = circum_loading_mid_dot = nothing
    end

    return OpenFASTData{TRadialInterpMethod,TTimeDerivMethod}(
        time, dtime_dtau, v, azimuth, omega, pitch,
        radii, radii_mid, dradii,
        cs_area, cs_area_mid,
        axial_loading, axial_loading_mid, axial_loading_mid_dot,
        circum_loading, circum_loading_mid, circum_loading_mid_dot)
end

function OpenFASTData{TRadialInterpMethod,TTimeDerivMethod}(
        time, dtime_dtau, v, azimuth, omega, pitch,
        radii, radii_mid, dradii,
        cs_area, cs_area_mid,
        axial_loading, axial_loading_mid, axial_loading_mid_dot,
        circum_loading, circum_loading_mid, circum_loading_mid_dot) where {TRadialInterpMethod,TTimeDerivMethod}
    return OpenFASTData{TRadialInterpMethod,TTimeDerivMethod,
               typeof(time),typeof(dtime_dtau),typeof(v),typeof(azimuth),typeof(omega),typeof(pitch),
               typeof(radii),typeof(radii_mid),typeof(dradii),
               typeof(cs_area),typeof(cs_area_mid),
               typeof(axial_loading),typeof(axial_loading_mid),typeof(axial_loading_mid_dot),
               typeof(circum_loading),typeof(circum_loading_mid),typeof(circum_loading_mid_dot)}(
                   time, dtime_dtau, v, azimuth, omega, pitch,
                   radii, radii_mid, dradii,
                   cs_area, cs_area_mid,
                   axial_loading, axial_loading_mid, axial_loading_mid_dot,
                   circum_loading, circum_loading_mid, circum_loading_mid_dot)
end

function _get_num_blades(fmt, column_names)
    # Match the format against the column names.
    ms = match.(fmt, column_names)

    # Remove any non-matches, for which `match` returns `nothing`.
    ms_only_matches = filter(x->!isnothing(x), ms)

    # Now, get all the sorted, unique blade indices, converting them from strings to Ints.
    blade_idxs = parse.(Int, unique(sort(getindex.(ms_only_matches, :blade))))

    # The blade indices appear to start at 1, so the highest blade index is the number of blades.
    num_blades = maximum(blade_idxs)

    # But check that the blade indices are what we assumed.
    @assert all(blade_idxs .== 1:num_blades)

    return num_blades, ms_only_matches
end

function _get_num_blades_num_radial(fmt, column_names)
    # Let's figure out how many blades and radial stations there are.
    # First, apply the loading regular expression to each column name:
    ms = match.(fmt, column_names)

    # Remove any non-matches, for which `match` returns `nothing`.
    ms_only_matches = filter(x->!isnothing(x), ms)

    # Now, get all the sorted, unique blade indices, converting them from strings to Ints.
    blade_idxs = parse.(Int, unique(sort(getindex.(ms_only_matches, :blade))))

    # The blade indices appear to start at 1, so the highest blade index is the number of blades.
    num_blades = maximum(blade_idxs)

    # But check that the blade indices are what we assumed.
    @assert all(blade_idxs .== 1:num_blades)

    # Now do the same thing for the radial indices.
    radial_idxs = parse.(Int, unique(sort(getindex.(ms_only_matches, :radial))))
    num_radial = maximum(radial_idxs)
    @assert all(radial_idxs .== 1:num_radial)

    return num_blades, num_radial, ms_only_matches
end

function interpolate_to_cell_centers!(data::OpenFASTData{FLOWLinearInterp})
    if (data.cs_area !== nothing) && (data.cs_area_mid !== nothing)
        data.cs_area_mid .= FLOWMath.linear.(Ref(data.radii), Ref(data.cs_area), data.radii_mid)
    end

    if (data.axial_loading !== nothing) && (data.axial_loading_mid !== nothing)
        for bidx in 1:size(data.axial_loading_mid, 3)
            for tidx in 1:size(data.axial_loading_mid, 1)
                data.axial_loading_mid[tidx, :, bidx] .= FLOWMath.linear.(Ref(data.radii), Ref(@view(data.axial_loading[tidx, :, bidx])), data.radii_mid)
            end
        end
    end

    if (data.circum_loading !== nothing) && (data.circum_loading_mid !== nothing)
        for bidx in 1:size(data.circum_loading_mid, 3)
            for tidx in 1:size(data.circum_loading_mid, 1)
                data.circum_loading_mid[tidx, :, bidx] .= FLOWMath.linear.(Ref(data.radii), Ref(@view(data.circum_loading[tidx, :, bidx])), data.radii_mid)
            end
        end
    end

    return nothing
end

function interpolate_to_cell_centers!(data::OpenFASTData{FLOWAkimaInterp})
    if (data.cs_area !== nothing) && (data.cs_area_mid !== nothing)
        spline_cs_area = FLOWMath.Akima(data.radii, data.cs_area)
        data.cs_area_mid .= spline_cs_area.(data.radii_mid)
    end

    if (data.axial_loading !== nothing) && (data.axial_loading_mid !== nothing)
        for bidx in 1:size(data.axial_loading_mid, 3)
            for tidx in 1:size(data.axial_loading_mid, 1)
                spline_axial = FLOWMath.Akima(data.radii, @view(data.axial_loading[tidx, :, bidx]))
                data.axial_loading_mid[tidx, :, bidx] .= spline_axial.(data.radii_mid)
            end
        end
    end

    if (data.circum_loading !== nothing) && (data.circum_loading_mid !== nothing)
        for bidx in 1:size(data.circum_loading_mid, 3)
            for tidx in 1:size(data.circum_loading_mid, 1)
                spline_circum = FLOWMath.Akima(data.radii, @view(data.circum_loading[tidx, :, bidx]))
                data.circum_loading_mid[tidx, :, bidx] .= spline_circum.(data.radii_mid)
            end
        end
    end

    return nothing
end

function calculate_loading_dot!(data::OpenFASTData{TRadialInterpMethod,NoTimeDerivMethod}) where {TRadialInterpMethod}
    if data.axial_loading_mid_dot !== nothing
        fill!(data.axial_loading_mid_dot, 0)
    end

    if data.circum_loading_mid_dot !== nothing
        fill!(data.circum_loading_mid_dot, 0)
    end

    return nothing
end

function _finite_diff_2nd_order!(df_dtau, f)
    # `f` is the input array, which we assume is of size `(num_times, num_radial, num_blades)`.
    # `df_dtau` is the derivative wrt `tau`, the non-dimensional time.

    @views begin
        # These stencils are in Tannehill, Anderson, Pletcher, "Computational Fluid Mechanics and Heat Transfer," 2nd edition, page 50.
        # First do the interior points.
        df_dtau[begin+1:end-1, :, :] .= 0.5 .* (f[begin+2:end, :, :] .- f[begin:end-2, :, :])

        # Then the left boundary.
        df_dtau[begin, :, :] .= 0.5 .* (-3 .* f[begin, :, :] .+ 4 .* f[begin+1, :, :] .- f[begin+2, :, :])

        # Then the right boundary.
        df_dtau[end, :, :] .= 0.5 .* (3 .* f[end, :, :] .- 4 .* f[end-1, :, :] .+ f[end-2, :, :])
    end

    return nothing
end

function calculate_loading_dot!(data::OpenFASTData{TRadialInterpMethod,SecondOrderFiniteDiff}) where {TRadialInterpMethod}
    # First get dt/dτ.
    _finite_diff_2nd_order!(data.dtime_dtau, data.time)

    if (data.axial_loading_mid !== nothing) && (data.axial_loading_mid_dot !== nothing)
        # Now get the derivatitve of the axial loading wrt tau, the non-dimensional time.
        _finite_diff_2nd_order!(data.axial_loading_mid_dot, data.axial_loading_mid)

        # Now get the derivative of the axial loading with respect to the dimensional time via `dfdt = df/dτ*dτ/dt = (df/dτ)/(dt/dτ)
        data.axial_loading_mid_dot ./= data.dtime_dtau
    end

    if (data.circum_loading_mid !== nothing) && (data.circum_loading_mid_dot !== nothing)
        # Now get the derivatitve of the circum loading wrt tau, the non-dimensional time.
        _finite_diff_2nd_order!(data.circum_loading_mid_dot, data.circum_loading_mid)

        # Now get the derivative of the circum loading with respect to the dimensional time via `dfdt = df/dτ*dτ/dt = (df/dτ)/(dt/dτ)
        data.circum_loading_mid_dot ./= data.dtime_dtau
    end

    return nothing
end

"""
    read_openfast_file(fname, radii, cs_area=nothing;
        header_keyword="Time",
        has_units_header=true,
        time_column_name="Time", 
        freestream_vel_column_name="Wind1VelX",
        azimuth_column_name="Azimuth",
        omega_column_name="RotSpeed",
        pitch_fmt=r"BlPitch(?<blade>[[:digit:]]+)",
        axial_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fxl",
        circum_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fyl",
        radial_interp_method=FLOWLinearInterp,
        time_deriv_method=SecondOrderFiniteDiff)

Read an OpenFAST output file and return a [`OpenFASTData`](@ref) object.

The `Azimuth` and `BlPitch` columns are assumed to be in degrees and will be converted to radians.
Likewise, the `RotSpeed` column is assumed to be in revolutions per minute and will be converted to radians per second.

# Arguments
* `fname`: name of the OpenFAST output file to read
* `radii`: `Vector` of blade radial coordinates
* `cs_area`: `Vector` of radial distribution of cross-sectional areas, or `nothing` to ignore
* `header_keyword="Time"`: string at the beginning of the header line (maybe always "Time"?)
* `has_units_header=true`: if true, assume the file has a line directly after the header line with the units of each column
* `time_column_name=header_keyword`: name of time column in file. Set to `nothing` to skip.
* `freestream_vel_column_name`: name of the freestream velocity column in the file. Set to `nothing` to skip.
* `azimuth_column_name`: name of the azimuth column in the file. Set to `nothing` to skip.
* `omega_column_name`: name of the omega (rotation rate) Set to `nothing` to skip.
* `pitch_fmt`: Format for finding all pitch columns in the file. Should be a regex with a capture group named `blade` for the blade index, or `nothing` to skip.
* `axial_loading_fmt`: Format for finding all axial loading columns in the file. Should be a regex with a captures groups named `blade` and `radial` for the blade and radial indices, or `nothing` to skip.
* `circum_loading_fmt`: Format for finding all radial loading columns in the file. Should be a regex with a captures groups named `blade` and `radial` for the blade and radial indices, or `nothing` to skip.
* `radial_interp_method`: `<:AbstractRadialInterpMethod` indicating method used to interpolate loading from blade element "interfaces" to midpoints.
* `time_deriv_method`: `<:AbstractTimeDerivMethod` indicating the method used to calculate the loading time derivatives.
* `average_freestream_vel=false`: Store possibily unsteady freestream velocity in the `OpenFASTData` object if `false`, store average value otherwise.
* `average_omega=false`: Store possibily unsteady omega (rotation rate) in the `OpenFASTData` object if `false`, store average value otherwise.
"""
function read_openfast_file(fname, radii, cs_area=nothing;
        header_keyword="Time",
        has_units_header=true,
        time_column_name=header_keyword, 
        freestream_vel_column_name="Wind1VelX",
        azimuth_column_name="Azimuth",
        omega_column_name="RotSpeed",
        pitch_fmt=r"BlPitch(?<blade>[[:digit:]]+)",
        axial_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fxl",
        circum_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fyl",
        radial_interp_method=FLOWLinearInterp,
        time_deriv_method=SecondOrderFiniteDiff,
        average_freestream_vel=true,
        average_omega=true)

    num_radial = length(radii)

    # Remove leading and trailing whitespace from header keyword.
    header_keyword_s = strip(header_keyword)

    # Find the first line that starts with `header_keyword_s`, which is where the header starts.
    idx_header = 0
    found_header = false
    for line in eachline(fname)
        idx_header += 1
        if startswith(strip(line), header_keyword_s)
            found_header = true
            break
        end
    end

    # If we didn't find a header, throw an error.
    if !found_header
        throw(ArgumentError("Unable to find header (line starting with \"$(header_keyword_s)\") in $(fname)"))
    end

    # Decide what to do with the units header.
    if has_units_header
        # If we have a units header, then we'll want to start reading right after it.
        skipto = idx_header + 2
    else
        # If there is no units header, then we don't need to skip anything and we'll just start reading directly after the header.
        skipto = idx_header + 1
    end

    # Read the file into a dataframe.
    df = CSV.read(fname, DataFrames.DataFrame; header=idx_header, skipto=skipto)

    # The number of times is equal to the number of rows in the dataframe.
    num_times = DataFrames.nrow(df)

    # This gives us all the column names in the dataframe.
    colnames = DataFrames.names(df)

    if time_column_name === nothing
        time = nothing
    else
        time = df[!, time_column_name]
    end

    if freestream_vel_column_name === nothing
        v = nothing
    else
        v_tmp = df[!, freestream_vel_column_name]
        if average_freestream_vel
            v = Fill(mean(v_tmp), length(v_tmp))
        else
            v = v_tmp
        end
    end

    if azimuth_column_name === nothing
        azimuth = nothing
    else
        azimuth = df[!, azimuth_column_name] .* (pi/180)
    end

    if omega_column_name === nothing
        omega = nothing
    else
        omega_tmp = df[!, omega_column_name] .* (2*pi/60)
        if average_omega
            omega = Fill(mean(omega_tmp), length(omega_tmp))
        else
            omega = omega_tmp
        end
    end

    if pitch_fmt === nothing
        pitch = nothing
    else
        # Get the number of blades and the pitch column names according to the pitch format.
        pitch_num_blades, pitch_matches = _get_num_blades(pitch_fmt, colnames)

        # Decide on an element type for the pitch, then read in the pitch.
        TF_pitch = promote_type(eltype.(getproperty.(Ref(df), getproperty.(pitch_matches, :match)))...)
        pitch = Array{TF_pitch, 2}(undef, num_times, pitch_num_blades)
        for m in pitch_matches
            b = parse(Int, m[:blade])
            pitch[:, b] .= df[!, m.match] .* (pi/180)
        end
    end

    if axial_loading_fmt === nothing
        axial_loading = nothing
    else
        # Get the number of blades and radial stations according to the axial loading format.
        axial_num_blades, axial_num_radial, axial_matches = _get_num_blades_num_radial(axial_loading_fmt, colnames)

        # Decide on an element type for the axial loading, then read in the axial loading.
        TF_axial = promote_type(eltype.(getproperty.(Ref(df), getproperty.(axial_matches, :match)))...)
        axial_loading = Array{TF_axial, 3}(undef, num_times, axial_num_radial, axial_num_blades)
        for m in axial_matches
            b = parse(Int, m[:blade])
            r = parse(Int, m[:radial])
            axial_loading[:, r, b] .= df[!, m.match]
        end
    end

    if circum_loading_fmt === nothing
        circum_loading = nothing
    else
        # Get the number of blades and radial stations according to the circumferential loading format.
        circum_num_blades, circum_num_radial, circum_matches = _get_num_blades_num_radial(circum_loading_fmt, colnames)

        # Decide on an element type for the circumferential loading.
        TF_circum = promote_type(eltype.(getproperty.(Ref(df), getproperty.(circum_matches, :match)))...)
        circum_loading = Array{TF_circum, 3}(undef, num_times, circum_num_radial, circum_num_blades)
        for m in circum_matches
            b = parse(Int, m[:blade])
            r = parse(Int, m[:radial])
            circum_loading[:, r, b] .= df[!, m.match]
        end
    end

    # Create the openfast data struct.
    data = OpenFASTData{radial_interp_method,time_deriv_method}(time, v, azimuth, omega, pitch, radii, cs_area, axial_loading, circum_loading)

    # Interpolate the loading to the cell centers.
    interpolate_to_cell_centers!(data)

    # Calculate the loading time derivatives.
    calculate_loading_dot!(data)

    return data
end

"""
    f1a_source_elements_openfast(data::OpenFASTData, rho0, c0, area_per_chord2::Vector, positive_x_rotation::Bool=true)

Construct and return an array of `CompactF1ASourceElement` objects from OpenFAST data.

# Arguments
- `data`: OpenFAST data object.
- `rho0`: Ambient air density (kg/m^3)
- `c0`: Ambient speed of sound (m/s)
- `positive_x_rotation`: rotate blade around the positive-x axis if `true`, negative-x axis otherwise
"""
function f1a_source_elements_openfast(data::OpenFASTData, rho0, c0, positive_x_rotation::Bool=true)

    # if length(area_per_chord2) != length(data.radii_mid)
    #     throw(ArgumentError("length of area_per_chord2 = $(length(area_per_chord2)) should be equal to length(data.radii_mid) = $(length(data.radii_mid))"))
    # end

    # OK, what's the coordinate system here?
    # Well, I guess I can decide.
    # For CCBlade, I've been assuming that the blade elements are translating in the positive x direction, rotating about the positive x axis if `positive_x_rotation` is `true`, negative x axis otherwise.
    # So for consistency let's say I do the same here.
    # Then that would mean the freestream velocity would be pointed in the negative x direction.
    # The velocity in the example OpenFAST file is always positive, so, I'll need to keep that in mind.
    # For the position of the blade elements, we'll initially be aligned with the y axis, and I'm assuming that the radial locations are all positive, so no sign switch necessary there.

    # But what to do about the loading?
    # What is the loading sign convention?
    # Well, for the axial loading, I would expect that the axial loading on the fluid would oppose the freestream velocity for a wind turbine.
    # So, if the blades are translating in the positive x direction, the axial loading on the fluid would also be in the positive x direction.
    # I can see that the axial loading in the openfast file is all positive, so no sign switch is necessary.
    # For the circumferential loading, for a wind turbine, I would expect the loading on the blade to be in the same direction as the blade rotation, so the loading on the fluid would be opposite the blade's rotation.
    # So, if the first blade is initially aligned with the y axis, rotating around the positive x axis, then the blade would initially be moving in the positive z direction, and so the loading on the fluid should be in the negative z direction.
    # And it looks like the circumferential loading in the OpenFAST file is negative, so no sign switch necessary.
    # But if the blade is rotating about the negative x axis, then the loading on the blade would be in the negative direction, and thus the loading on the fluid would be in the positive direction, and so we'd need to switch the sign.

    # Now, let's get some transformation stuff going.
    # We're going to be translating in the positive x direction.
    # Oh, but what do we do about the fact that the axial velocity isn't necessarily constant?
    # Well, what I need is the position and velocity in the axial direction.
    # It would also be nice to get the acceleration and jerk, too.
    # OK, I have the velocity, but the position is the tricky part.
    # But I should be able to just integrate it, I suppose.
    # So let's do that.
    # We'll assume the position at the first time is at the origin.
    # Now we need to integrate the velocity.
    # (This assumes the hub starts at the origin.
    # If that's not the case, just add `x0` to it, where `x0` is the axial position at the first time level.)
    x = FLOWMath.cumtrapz(data.time, data.v)

    # So now we have the axial position and axial velocity as a function of data.time.
    # Now we should be able to create transformations that take that into account.
    # There will be one for each data.time level.
    # For each entry in `data.time`, i.e. for `data.time[i]`, it will start at `[x[i], 0, 0]` and have velocity `[v[i], 0, 0]`.
    # This won't take into account the effect the possibily non-constant velocity has on the acceleration or jerk.
    # const_vel_trans = ConstantVelocityTransformation.(data.time, SVector{3,typeof(x)}.(x, 0, 0), SVector{3,typeof(x)}.(data.v, 0, 0))

    # Next, we need to figure out the rotation stuff.
    # We actually have both a time history of omega and azimuthal angles.
    # So we could create a bunch of rotational transformations that have the correct azimuth offset and omega.
    # But it wouldn't take into account the effect the time derivative of omega has on the things we care about: the derivatives of position and loading.
    # I would need to do some work on that.
    # Also, I'm going to assume that omega and azimuth are always positive, and so I'll switch their signs if we are rotating about the negative x axis.
    # rot_trans = SteadyRotXTransformation.(data.time, data.omega.*ifelse(positive_x_rotation, 1, -1), data.azimuth.*ifelse(positive_x_rotation, 1, -1))

    # So, want everything to have shape (num_times, num_radial_mid, num_blades).
    radii_mid_rs = reshape(data.radii_mid, 1, :)
    dradii_rs = reshape(data.dradii, 1, :)

    num_blades = data.num_blades
    thetas = 2*pi/num_blades.*(0:(num_blades-1)) .* ifelse(positive_x_rotation, 1, -1)
    thetas_rs = reshape(thetas, 1, 1, :)
    # cs_area_rs = reshape(cs_area, 1, :)
    cs_area_rs = reshape(data.cs_area_mid, 1, :)

    # All the loading arrays are in the right shape, and no sign switch is necessary.
    fn = data.axial_loading_mid
    fndot = data.axial_loading_mid_dot
    # The OpenFAST file doesn't have any data for the loading in the radial direction (which should be quite small of course).
    fr = zero(eltype(fn))
    frdot = zero(eltype(fn))
    fc = data.circum_loading_mid
    fcdot = data.circum_loading_mid_dot

    # Now we should be able to do all this.
    # This is a bit too cute, but I can't help myself.
    ses = (CompactF1ASourceElement.(rho0, c0, radii_mid_rs, thetas_rs, dradii_rs, cs_area_rs, fn, fndot, fr, frdot, fc.*ifelse(positive_x_rotation, 1, -1), fcdot.*ifelse(positive_x_rotation, 1, -1), data.time) .|>
           compose.(data.time,
                    ConstantVelocityTransformation.(data.time, SVector{3,eltype(x)}.(x, 0, 0), SVector{3,eltype(data.v)}.(data.v, 0, 0)),
                    SteadyRotXTransformation.(data.time, data.omega.*ifelse(positive_x_rotation, 1, -1), data.azimuth.*ifelse(positive_x_rotation, 1, -1))
                   )
          )

    return ses
end
