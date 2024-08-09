"""
Struct for holding data from an OpenFAST (AeroDyn?) output file.

# Fields:
* `time`: vector of simulation times with size `(num_times,)`
* `v`: vector of freestream velocity time history with size `(num_times,)`
* `azimuth`: vector of azimuth angle time history with size `(num_times,)`
* `omega`: vector of rotation rate time history with size `(num_times,)`
* `pitch`: array of pitch angle time history with size `(num_times, num_blades)`
* `axial_loading`: array of axial loading time history with size `(num_times, num_radial, num_blades)`
* `circum_loading`: array of circumferential loading time history with size `(num_times, num_radial, num_blades)`
"""
struct OpenFASTData{TTime,TV,TAzimuth,TOmega,TPitch,TAxialLoading,TTangentialLoading}
    time::TTime
    v::TV
    azimuth::TAzimuth
    omega::TOmega
    pitch::TPitch
    axial_loading::TAxialLoading
    circum_loading::TTangentialLoading

    function OpenFASTData(time, v, azimuth, omega, pitch, axial_loading, circum_loading)
        # Figure out what num_times is.
        if time !== nothing
            num_times = length(time)
        elseif v !== nothing
            num_times = length(v)
        elseif azimuth !== nothing
            num_times = length(azimuth)
        elseif omega !== nothing
            num_times = length(omega)
        elseif pitch !== nothing
            num_times = size(pitch, 1)
        elseif axial_loading !== nothing
            num_times = size(axial_loading, 1)
        elseif circum_loading !== nothing
            num_times = size(circum_loading, 1)
        end

        # Figure out what num_blades is.
        if pitch !== nothing
            num_blades = size(pitch, 2)
        elseif axial_loading !== nothing
            num_blades = size(axial_loading, 3)
        elseif circum_loading !== nothing
            num_blades = size(circum_loading, 3)
        end

        # Figure out what num_radial is.
        if axial_loading !== nothing
            num_radial = size(axial_loading, 2)
        elseif circum_loading !== nothing
            num_radial = size(circum_loading, 2)
        end

        if time !== nothing
            # I don't think this can ever happen since the check for time !== nothing at the beginning is first, but whatever.
            size(time) == (num_times,) || raise(ArgumentError("size(time) = $(size(time)) inconsistent with other inputs"))
        end

        if v !== nothing
            size(v) == (num_times,) || raise(ArgumentError("size(v) = $(size(v)) does not match size(time) = $(size(time))"))
        end

        if azimuth !== nothing
            size(azimuth) == (num_times,) || raise(ArgumentError("size(azimuth) = $(size(azimuth)) does not match size(time) = $(size(time))"))
        end

        if omega !== nothing
            size(omega) == (num_times,) || raise(ArgumentError("size(omega) = $(size(omega)) does not match size(time) = $(size(time))"))
        end

        if pitch !== nothing
            size(pitch) == (num_times, num_blades) || raise(ArgumentError("size(pitch) = $(size(pitch)) not consistent with size(time) = $(size(time)) and/or size(axial_loading) = $(size(axial_loading))"))
        end

        if axial_loading !== nothing
            size(axial_loading) == (num_times, num_radial, num_blades) || raise(ArgumentError("size(axial_loading) = $(size(axial_loading)) not consistent with size(time) = $(size(time))"))
        end

        if circum_loading !== nothing
            size(circum_loading) == (num_times, num_radial, num_blades) || raise(ArgumentError("size(circum_loading) = $(size(circum_loading)) not consistent with size(time) = $(size(time)) and/or size(axial_loading) = $(size(axial_loading))"))
        end

        return new{typeof(time),typeof(v),typeof(azimuth),typeof(omega),typeof(pitch),typeof(axial_loading),typeof(circum_loading)}(time, v, azimuth, omega, pitch, axial_loading, circum_loading)
    end
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

"""
    read_openfast_file(fname;
        header_keyword="Time",
        has_units_header=true,
        time_column_name="Time", 
        freestream_vel_column_name="Wind1VelX",
        azimuth_column_name="Azimuth",
        omega_column_name="RotSpeed",
        pitch_fmt=r"BlPitch(?<blade>[[:digit:]]+)",
        axial_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fxl",
        circum_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fyl")

Read an OpenFAST output file and return a [`OpenFASTData`](@ref) object.

# Arguments
* `fname`: name of the OpenFAST output file to read
* `header_keyword="Time"`: string at the beginning of the header line (maybe always "Time"?)
* `has_units_header=true`: if true, assume the file has a line directly after the header line with the units of each column
* `time_column_name=header_keyword`: name of time column in file. Set to `nothing` to skip.
* `freestream_vel_column_name`: name of the freestream velocity column in the file. Set to `nothing` to skip.
* `azimuth_column_name`: name of the azimuth column in the file. Set to `nothing` to skip.
* `omega_column_name`: name of the omega column in the file. Set to `nothing` to skip.
* `pitch_fmt`: Format for finding all pitch columns in the file. Should be a regex with a capture group named `blade` for the blade index, or `nothing` to skip.
* `axial_loading_fmt`: Format for finding all axial loading columns in the file. Should be a regex with a captures groups named `blade` and `radial` for the blade and radial indices, or `nothing` to skip.
* `circum_loading_fmt`: Format for finding all radial loading columns in the file. Should be a regex with a captures groups named `blade` and `radial` for the blade and radial indices, or `nothing` to skip.
"""
function read_openfast_file(fname;
        header_keyword="Time",
        has_units_header=true,
        # time_column_name="Time", 
        time_column_name=header_keyword, 
        freestream_vel_column_name="Wind1VelX",
        azimuth_column_name="Azimuth",
        omega_column_name="RotSpeed",
        pitch_fmt=r"BlPitch(?<blade>[[:digit:]]+)",
        axial_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fxl",
        circum_loading_fmt=r"AB(?<blade>[[:digit:]]+)N(?<radial>[[:digit:]]+)Fyl")

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
        v = df[!, freestream_vel_column_name]
    end

    if azimuth_column_name === nothing
        azimuth = nothing
    else
        azimuth = df[!, azimuth_column_name]
    end

    if omega_column_name === nothing
        omega = nothing
    else
        omega = df[!, omega_column_name]
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
            pitch[:, b] .= df[!, m.match]
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
        # Make sure we get the same thing if we use the circumferential loading format.
        circum_num_blades, circum_num_radial, circum_matches = _get_num_blades_num_radial(circum_loading_fmt, colnames)

        # Decide on an element type for the axial loading.
        TF_circum = promote_type(eltype.(getproperty.(Ref(df), getproperty.(circum_matches, :match)))...)
        circum_loading = Array{TF_circum, 3}(undef, num_times, circum_num_radial, circum_num_blades)
        for m in circum_matches
            b = parse(Int, m[:blade])
            r = parse(Int, m[:radial])
            circum_loading[:, r, b] .= df[!, m.match]
        end
    end

    return OpenFASTData(time, v, azimuth, omega, pitch, axial_loading, circum_loading)
end
