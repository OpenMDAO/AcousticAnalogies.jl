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
    lines = [MeshCell(PolyData.Lines(), line) for line in eachcol(line_ids)]

    vtkfile = vtk_grid(name, points, lines)

    _write_data_to_vtk!(vtkfile, ses)

    return vtkfile
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:CompactSourceElement})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    vtkfile["Length", VTKCellData()] = mapview(:Δr, ses)
    vtkfile["CSArea", VTKCellData()] = mapview(:Λ, ses)
    vtkfile["Position", VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["Acceleration", VTKCellData()] = hcat(mapview(:y2dot, ses)...)
    vtkfile["Jerk", VTKCellData()] = hcat(mapview(:y3dot, ses)...)
    vtkfile["Loading", VTKCellData()] = hcat(mapview(:f0dot, ses)...)
    vtkfile["LoadingDot", VTKCellData()] = hcat(mapview(:f1dot, ses)...)
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:Union{TBLTESourceElement,LBLVSSourceElement,TipVortexSourceElement}})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    # vtkfile["Length", VTKCellData()] = SingleFieldStructArray(ses, Val{:Δr})
    # vtkfile["Chord", VTKCellData()] = SingleFieldStructArray(ses, Val{:chord})
    # vtkfile["Position", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y0dot})...)
    # vtkfile["Velocity", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y1dot})...)
    # vtkfile["FluidVelocity", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y1dot_fluid})...)
    # vtkfile["SpanUnitVector", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:span_uvec})...)
    # vtkfile["ChordUnitVector", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:chord_uvec})...)

    vtkfile["Length", VTKCellData()] = mapview(:Δr, ses)
    vtkfile["Chord", VTKCellData()] = mapview(:chord, ses)
    vtkfile["Position", VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["FluidVelocity", VTKCellData()] = hcat(mapview(:y1dot_fluid, ses)...)
    vtkfile["SpanUnitVector", VTKCellData()] = hcat(mapview(:span_uvec, ses)...)
    vtkfile["ChordUnitVector", VTKCellData()] = hcat(mapview(:chord_uvec, ses)...)
end

function _write_data_to_vtk!(vtkfile, ses::AbstractArray{<:Union{TEBVSSourceElement,CombinedNoTipBroadbandSourceElement,CombinedWithTipBroadbandSourceElement}})
    # Now need to add the cell data. I would have expected to have to flatten
    # these arrays, but apparently that's not necessary.
    # vtkfile["Length", VTKCellData()] = SingleFieldStructArray(ses, Val{:Δr})
    # vtkfile["Chord", VTKCellData()] = SingleFieldStructArray(ses, Val{:chord})
    # vtkfile["TrailingEdgeThickness", VTKCellData()] = SingleFieldStructArray(ses, Val{:h})
    # vtkfile["TrailingEdgeAngle", VTKCellData()] = SingleFieldStructArray(ses, Val{:Psi})
    # vtkfile["Position", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y0dot})...)
    # vtkfile["Velocity", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y1dot})...)
    # vtkfile["FluidVelocity", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:y1dot_fluid})...)
    # vtkfile["SpanUnitVector", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:span_uvec})...)
    # vtkfile["ChordUnitVector", VTKCellData()] = hcat(SingleFieldStructArray(ses, Val{:chord_uvec})...)

    vtkfile["Length", VTKCellData()] = mapview(:Δr, ses)
    vtkfile["Chord", VTKCellData()] = mapview(:chord, ses)
    vtkfile["TrailingEdgeThickness", VTKCellData()] = mapview(:h, ses)
    vtkfile["TrailingEdgeAngle", VTKCellData()] = mapview(:Psi, ses)
    vtkfile["Position", VTKCellData()] = hcat(mapview(:y0dot, ses)...)
    vtkfile["Velocity", VTKCellData()] = hcat(mapview(:y1dot, ses)...)
    vtkfile["FluidVelocity", VTKCellData()] = hcat(mapview(:y1dot_fluid, ses)...)
    vtkfile["SpanUnitVector", VTKCellData()] = hcat(mapview(:span_uvec, ses)...)
    vtkfile["ChordUnitVector", VTKCellData()] = hcat(mapview(:chord_uvec, ses)...)
end

function to_paraview_collection!(pvd, name, ses::AbstractArray{<:AbstractCompactSourceElement}, time_axis=1)
    # These are the time indices, i.e., the ones we actually want to iterate
    # over.
    indices = axes(ses, time_axis)

    # idx is all the indicies, including the time one that we want to iterate
    # over, and the others that we want to grab all at once. They're all colons
    # right now, but the time one will be replaced with each index in `indices`
    # in the `for i in indices` loop below.
    idx = Any[Colon() for ind in axes(ses)]

    for i in indices
        idx[time_axis] = i
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

"""
    to_paraview_collection(name::AbstractString, ses::NTuple{N, AbstractArray{<:AbstractCompactSourceElement}}, time_axes::NTuple{N, Int64}=ntuple(i->1, N))


Construct and write out a ParaView collection data file (`.pvd`) object for a tuple of arrays of `CompactSourceElement`s with name `name.pvd` (i.e., the `name` argument should not contain a file extension).

`time_axes` is a tuple of time axis indices, one for each entry of `ses`, indicating the axis over which the source time for the source elements in `ses` vary.
One VTK PolyData (`.vtp`) file will be written for each valid index along `time_axis` for each array in the `ses` tuple.

Returns a list of filenames written out by [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).
"""
function to_paraview_collection(name, ses::NTuple{N, AbstractArray{<:AbstractCompactSourceElement}}, time_axes::NTuple{N, Int64}=ntuple(i->1, N)) where {N}

    outfiles = paraview_collection(name) do pvd
        for (i, (sesi, time_axis)) in enumeratet(zip(ses, time_axes))
            namei = format("{}-{:03d}-", name, i)
            to_paraview_collection!(pvd, namei, sesi, time_axis)
        end
    end

    return outfiles
end

"""
    to_paraview_collection(name::AbstractString, ses::AbstractArray{<:AbstractCompactSourceElement}, time_axis::Integer=1)


Construct and write out a ParaView collection data file (`.pvd`) object for an array of `CompactSourceElement` with name `name.pvd` (i.e., the `name` argument should not contain a file extension).

`time_axis` indicates the time_axis of `ses` over which the source time for the source
elements in `ses` vary. One VTK PolyData (`.vtp`) file will be written for each
valid index along `time_axis`.

Returns a list of filenames written out by [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl).
"""
function to_paraview_collection(name, ses::AbstractArray{<:AbstractCompactSourceElement}, time_axis=1)

    outfiles = paraview_collection(name) do pvd
        to_paraview_collection!(pvd, name, ses, time_axis)
    end

    return outfiles
end
