"""
    endpoints(se::AbstractCompactSourceElement)

Return the Tuple containing the endpoint locations of the compact source element `se`.
"""
function endpoints(se::AbstractCompactSourceElement)
    p1 = se.y0dot - 0.5*se.Δr*orientation(se)
    p2 = se.y0dot + 0.5*se.Δr*orientation(se)
    # Don't like the idea of having a billion SVectors, so convert them to plain
    # Vectors. Maybe not necessary.
    return (Vector(p1), Vector(p2))
end

"""
    to_vtp(name::AbstractString, ses::AbstractArray{<:AbstractCompactSourceElement})

Construct and return a VTK polygonal (.vtp) data file object for an array of `AbstractCompactSourceElement` with name `name.vtp` (i.e., the `name` argument should not contain a file extension).
"""
function to_vtp(name, ses::AbstractArray{<:AbstractCompactSourceElement})
    # This will be an array of the same size as ses, with each entry a
    # length-two tuple of the endpoints.
    points_all = AcousticAnalogies.endpoints.(ses)
    # This is an array of shape (2, size(ses)...) that has the global ID of each
    # point. Need to use `collect` to create an actual Array since we're going
    # to modify its elements.
    line_ids = reshape(1:(2*length(ses)), 2, size(ses)...) |> collect

    # Now I need to loop over each point...
    points = Vector{Vector{Float64}}()
    for I in CartesianIndices(points_all)
        for (j, p_current) in enumerate(points_all[I])
            # Is the current point in points?
            k = findfirst(p -> all(p .≈ p_current), points)
            if k === nothing
                # p_current is not in points, so add it.
                push!(points, p_current)
                k = length(points)
            end
            line_ids[j, I] = k
        end
    end

    # This converts the points Vector{Vector} into a Matrix. Should be
    # size (num_dims, num_unique_points) where `num_dims` should be 3 (number of
    # spatial dimensions), and num_unique_points is the number of unique points.
    points = hcat(points...)

    # `line_ids` is an Array of shape (2, size(ses)...) Need to reshape that to
    # (2, length(ses)).
    line_ids = reshape(line_ids, 2, length(ses))
    # Now create the array of VTK lines.
    lines = [WriteVTK.MeshCell(WriteVTK.PolyData.Lines(), line) for line in eachcol(line_ids)]

    vtkfile = WriteVTK.vtk_grid(name, points, lines)

    _write_data_to_vtk!(vtkfile, ses)

    return vtkfile
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:CompactF1ASourceElement})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    vtkfile["Length", WriteVTK.VTKCellData()] = mapview(:Δr, ses)
    vtkfile["CSArea", WriteVTK.VTKCellData()] = mapview(:Λ, ses)
    vtkfile["Position", WriteVTK.VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", WriteVTK.VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["Acceleration", WriteVTK.VTKCellData()] = hcat(mapview(:y2dot, ses)...)
    vtkfile["Jerk", WriteVTK.VTKCellData()] = hcat(mapview(:y3dot, ses)...)
    vtkfile["Loading", WriteVTK.VTKCellData()] = hcat(mapview(:f0dot, ses)...)
    vtkfile["LoadingDot", WriteVTK.VTKCellData()] = hcat(mapview(:f1dot, ses)...)
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:Union{TBLTESourceElement,LBLVSSourceElement,TipVortexSourceElement}})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    vtkfile["Length", WriteVTK.VTKCellData()] = mapview(:Δr, ses)
    vtkfile["Chord", WriteVTK.VTKCellData()] = mapview(:chord, ses)
    vtkfile["Position", WriteVTK.VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", WriteVTK.VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["FluidVelocity", WriteVTK.VTKCellData()] = hcat(mapview(:y1dot_fluid, ses)...)
    vtkfile["SpanUnitVector", WriteVTK.VTKCellData()] = hcat(mapview(:span_uvec, ses)...)
    vtkfile["ChordUnitVector", WriteVTK.VTKCellData()] = hcat(mapview(:chord_uvec, ses)...)
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:Union{TEBVSSourceElement,CombinedNoTipBroadbandSourceElement,CombinedWithTipBroadbandSourceElement}})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    vtkfile["Length", WriteVTK.VTKCellData()] = mapview(:Δr, ses)
    vtkfile["Chord", WriteVTK.VTKCellData()] = mapview(:chord, ses)
    vtkfile["TrailingEdgeThickness", WriteVTK.VTKCellData()] = mapview(:h, ses)
    vtkfile["TrailingEdgeAngle", WriteVTK.VTKCellData()] = mapview(:Psi, ses)
    vtkfile["Position", WriteVTK.VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", WriteVTK.VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["FluidVelocity", WriteVTK.VTKCellData()] = hcat(mapview(:y1dot_fluid, ses)...)
    vtkfile["SpanUnitVector", WriteVTK.VTKCellData()] = hcat(mapview(:span_uvec, ses)...)
    vtkfile["ChordUnitVector", WriteVTK.VTKCellData()] = hcat(mapview(:chord_uvec, ses)...)
end

function to_vtu(name, obs::AbstractAcousticObserver, t; sphere_radius=1.0, n=10)
    # Get the current position of the observer.
    x = obs(t)

    # Create a sphere at the observer position with the specified radius.
    s = Meshes.Sphere(tuple(x...), sphere_radius)

    # Discretize the sphere.
    mesh = Meshes.simplexify(Meshes.discretize(s, Meshes.RegularDiscretization(n)))

    # This gets a fancy Unitful vector of Points that represent the verticies of the discretized sphere.
    verts = Meshes.vertices(mesh)

    # This gets the connectivity of the discretized sphere.
    connec = Meshes.elements(Meshes.topology(mesh))

    # This converts the fancy `verts` into a plain old Matrix{Float64}.
    points = stack(p -> Meshes.ustrip.(Meshes.to(p)), verts)

    # This creates VTK cells that we can eventually write out, converting the Meshes.jl connectivity into VTK MeshCells.
    cells = [WriteVTK.MeshCell(WriteVTK.VTKCellTypes.VTK_TRIANGLE, Meshes.indices(c)) for c in connec]

    # Now we can finally create a VTK file.
    vtkfile = WriteVTK.vtk_grid(name, points, cells)

    return vtkfile
end

r"""
to_paraview_collection(name::AbstractString, ses::NTuple{N, AbstractArray{<:AbstractCompactSourceElement}}; time_axes::NTuple{N, Int64}=ntuple(i->1, N), block_names=["block$(b)" for b in 1:N], observers=(), observer_names=nothing, observer_radii=nothing)


Construct and write out a ParaView collection data file (`.pvd`) object for a tuple of arrays of `CompactF1ASourceElement`s with name `name.pvd` (i.e., the `name` argument should not contain a file extension).

`time_axes` is a tuple of time axis indices, one for each entry of `ses`, indicating the axis over which the source time for the source elements in `ses` vary.
`block_names` is a `Vector` of strings, one for each entry in `ses`, that will be used to name each `ses` entry in the multiblock VTK file.

`observers` should be an iterable of `AbstractAcousticObserver` objects that will also be written out to the paraview file as discretized spheres.
`observer_names` can either be an iterable of Strings that will be used to name each VTK observer file, or `nothing`, in which case each observer will be named `observer<int>`.
`observer_radii` can either be an iterable of Float64 representing the radius of each observer sphere, or `nothing`, in which case a suitible radius will be calculated.

One VTK PolyData (`.vtp`) file will be written for each valid index along `time_axis` for each array in the `ses` tuple.

Returns a list of filenames written out by [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).
"""
function to_paraview_collection(name, ses::NTuple{N, AbstractArray{<:AbstractCompactSourceElement}}; time_axes::NTuple{N, Int64}=ntuple(i->1, N), block_names=["block$(b)" for b in 1:N], observers=(), observer_names=nothing, observer_radii=nothing) where {N}

    # Check that the size of each array in the `ses` is the the same along the time axis.
    # This will get the length of each array in `ses` along the time axis.
    time_axis_lengths = size.(ses, time_axes)
    # Now check if they're the same.
    if !all(time_axis_lengths .== first(time_axis_lengths))
        throw(ArgumentError("length of each array in ses tuple must be the same along the time axis. Length of each: $(time_axis_lengths)"))
    end
    time_axis_len = first(time_axis_lengths)

    # Also check that the block names are unique.
    length(unique(block_names)) == N || throw(ArgumentError("entries in block_names = $(block_names) must be unique"))

    if length(observers) > 0
        if observer_names === nothing
            obs_names = ["observer$(i)" for i in 1:length(observers)]
        else
            # Check that the number of observers matches the number of observer names.
            length(observer_names) == length(observers) || throw(ArgumentError("length of observer_names does not match length of observers"))
            # Check that the observer names are unique.
            length(unique(observer_names)) == length(observer_names) || throw(ArgumentError("entries in observer_names = $(observer_names) must be unique"))

            # Also should check that the observer names aren't the same as the source element block names.
            all_names = vcat(block_names, observer_names)
            length(unique(all_names)) == length(all_names) || throw(ArgumentError("observer_names = $(observer_names) shares entries with block_names = $(block_names)"))

            obs_names = observer_names
        end

        if observer_radii === nothing
            # Check out this sorcery.
            # (xmin, xmax) = mapreduce(x->extrema(getindex.(getproperty.(x, :y0dot), 1)), (a, b)->(min(a[1], b[1]), max(a[2], b[2])), ses)
            # (ymin, ymax) = mapreduce(x->extrema(getindex.(getproperty.(x, :y0dot), 2)), (a, b)->(min(a[1], b[1]), max(a[2], b[2])), ses)
            # (zmin, zmax) = mapreduce(x->extrema(getindex.(getproperty.(x, :y0dot), 3)), (a, b)->(min(a[1], b[1]), max(a[2], b[2])), ses)
            # But could I do it in one line?
            # That would mean I'd need to have a function that could take in a vector, then compare it to another vector.
            # ((xmin, ymin, zmin), (xmax, ymax, zmax)) = mapreduce(
            #     (a, b)->((min(a[1][1], b[1][1]),
            #               min(a[1][2], b[1][2]),
            #               min(a[1][3], b[1][3])),
            #              (max(a[2][1], b[2][1]),
            #               max(a[2][2], b[2][2]),
            #               max(a[2][3], b[2][3]))), ses) do x
            #     xmin, xmax = extrema(getindex.(getproperty.(x, :y0dot), 1))
            #     ymin, ymax = extrema(getindex.(getproperty.(x, :y0dot), 2))
            #     zmin, zmax = extrema(getindex.(getproperty.(x, :y0dot), 3))
            #     return ((xmin, ymin, zmin), (xmax, ymax, zmax))
            # end
            # Still kind of annoying that I have to iterate over each source element array three times.
            # How to avoid that?
            # It would come down to not using extrema.
            # I'd need a version of extrema that iterates over each vector, keeping track of the minimum and maximum of each element.
            ((xmin, ymin, zmin), (xmax, ymax, zmax)) = mapreduce(
                (a, b)->((min(a[1][1], b[1][1]),
                          min(a[1][2], b[1][2]),
                          min(a[1][3], b[1][3])),
                         (max(a[2][1], b[2][1]),
                          max(a[2][2], b[2][2]),
                          max(a[2][3], b[2][3]))), ses) do X
                
                mins = (Inf, Inf, Inf)
                maxs = (-Inf, -Inf, -Inf)
                for x in X
                    y0dot = x.y0dot
                    mins = (min(y0dot[1], mins[1]), min(y0dot[2], mins[2]), min(y0dot[3], mins[3]))
                    maxs = (max(y0dot[1], maxs[1]), max(y0dot[2], maxs[2]), max(y0dot[3], maxs[3]))
                end
                return mins, maxs
            end

            # Now find the diagonal of the boundary box.
            r = sqrt((xmax - xmin)^2 + (ymax - ymin)^2 + (zmax - zmin)^2)

            # Set the observer radii to be a fraction of that boundary box diagonal
            obs_radii = fill(0.02*r, length(observers))
        else
            # Check that the length of observer_radii matches observers.
            length(observer_radii) == length(observers) || throw(ArgumentError("length of observer_radii does not match length of observers"))
            obs_radii = observer_radii
        end
    else
        obs_names = nothing
        obs_radii = nothing
    end

    outfiles = WriteVTK.paraview_collection(name) do pvd
        for tidx in 1:time_axis_len
            # We want to check that the source times for all elements we write out for this time index `tidx` is the same.
            # So get the time for the first element in the first array, and compare later.
            ses_first = first(ses)
            time_axis_first = first(time_axes)
            axes_first = axes(ses_first)
            idx_first = [d == time_axis_first ? axes_first[d][tidx] : axes_first[d][begin] for d in 1:ndims(ses_first)]
            t = source_time(ses_first[idx_first...])

            # We will create a multiblock VTK file for each time step, with one block for each array in `ses`.
            name_mb = format("{}-{:08d}", name, tidx)
            WriteVTK.vtk_multiblock(name_mb) do vtm
                for (i, (sesi, time_axis)) in enumerate(zip(ses, time_axes))
                    # This will give me each valid index allong the time axis for `sesi`.
                    time_indices = axes(sesi, time_axis)

                    # This will give me a `:` for each dimension of `sesi` except for the dimension corresponding to `time_axis`.
                    idx = [d == time_axis ? time_indices[tidx] : Colon() for d in 1:ndims(sesi)]

                    # Now we can grab all the source elements we want.
                    sesiv = @view sesi[idx...]

                    # Now check if all the source times are what we expect.
                    all(source_time.(sesiv) .≈ t) || thow(ArgumentError("ses[$i][$(idx)] does not have all source elements with source time $(t)"))

                    # Now create a VTK file for the source elements we've selected.
                    namei = format("{}-{}-{:08d}", name, block_names[i], tidx)
                    vtkfile = to_vtp(namei, sesiv)

                    # And add it to the multiblock file.
                    WriteVTK.multiblock_add_block(vtm, vtkfile)
                end

                for (obs, obs_name, obs_radius) in zip(observers, obs_names, obs_radii)
                    # Now also need to write out the observers.
                    namei = format("{}-{}-{:08d}", name, obs_name, tidx)
                    vtkfile = to_vtu(namei, obs, t; sphere_radius=obs_radius)

                    WriteVTK.multiblock_add_block(vtm, vtkfile)
                end

                # Now add the the multiblock file to the paraview collection.
                pvd[t] = vtm
            end
        end
    end

    return outfiles
end

"""
    to_paraview_collection(name::AbstractString, ses::AbstractArray{<:AbstractCompactSourceElement}; time_axis::Integer=1)


Construct and write out a ParaView collection data file (`.pvd`) object for an array of `AbstractCompactSourceElement`s with name `name.pvd` (i.e., the `name` argument should not contain a file extension).

`time_axis` indicates the time_axis of `ses` over which the source time for the source
elements in `ses` vary. One VTK PolyData (`.vtp`) file will be written for each
valid index along `time_axis`.

Returns a list of filenames written out by [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).
"""
function to_paraview_collection(name, ses::AbstractArray{<:AbstractCompactSourceElement}; time_axis=1)

    outfiles = WriteVTK.paraview_collection(name) do pvd
        # These are the time indices, i.e., the ones we actually want to iterate
        # over.
        time_indices = axes(ses, time_axis)

        # idx is all the indicies, including the time one that we want to iterate
        # over, and the others that we want to grab all at once. They're all colons
        # right now, but the time one will be replaced with each index in `indices`
        # in the `for i in indices` loop below.
        idx = [d == time_axis ? first(time_indices) : Colon() for d in 1:ndims(ses)]

        for (i, it) in enumerate(time_indices)
            idx[time_axis] = it
            namei = format("{}{:08d}", name, i)
            sesi = @view ses[idx...]
            vtkfile = to_vtp(namei, sesi)

            # Assume that all the times in sesi are the same. They should be if the
            # user specified the correct time_axis. I suppose I could check that and
            # issue a warning.
            t = first(sesi).τ
            pvd[t] = vtkfile
        end
    end

    return outfiles
end
