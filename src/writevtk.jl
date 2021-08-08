"""
    endpoints(se::CompactSourceElement)

Return the Tuple containing the endpoint locations of the compact source element `se`.
"""
function endpoints(se::CompactSourceElement)
    p1 = se.y0dot - 0.5*se.Δr*se.u
    p2 = se.y0dot + 0.5*se.Δr*se.u
    # Don't like the idea of having a billion SVectors, so convert them to plain
    # Vectors. Maybe not necessary.
    return (Vector(p1), Vector(p2))
end

function to_vtp(name, ses::AbstractArray{<:CompactSourceElement})
    # First need to get the nodes we'll use for the VTK file. This flattens a
    # Vector{Tuple{Vector, Vector}} into a Vector{Vector}.
    points = [(AcousticAnalogies.endpoints.(ses)...)...]
    # And this converts the Vector{Vector} into a Matrix.
    points = hcat(points...)

    lines = [MeshCell(PolyData.Lines(), (i, i+1)) for i in 1:2:(2*length(ses))]
    vtkfile = vtk_grid(name, points, lines)

    # Now need to add the cell data.
    vtkfile["Length", VTKCellData()] = SingleFieldStructArray(ses, :Δr)
    vtkfile["CSArea", VTKCellData()] = SingleFieldStructArray(ses, :Λ)
    vtkfile["Position", VTKCellData()] = hcat(SingleFieldStructArray(ses, :y0dot)...)
    vtkfile["Velocity", VTKCellData()] = hcat(SingleFieldStructArray(ses, :y1dot)...)
    vtkfile["Acceleration", VTKCellData()] = hcat(SingleFieldStructArray(ses, :y2dot)...)
    vtkfile["Jerk", VTKCellData()] = hcat(SingleFieldStructArray(ses, :y3dot)...)
    vtkfile["Loading", VTKCellData()] = hcat(SingleFieldStructArray(ses, :f0dot)...)
    vtkfile["LoadingDot", VTKCellData()] = hcat(SingleFieldStructArray(ses, :f1dot)...)
    return vtkfile
end

function to_paraview_collection(name, ses, axis=1)
    pvd = paraview_collection(name)

    # These are the time indices, i.e., the ones we actually want to iterate
    # over.
    indices = axes(ses, axis)

    # idx is all the indicies, including the time one that we want to iterate
    # over, and the others that we want to grab all at once. They're all colons
    # right now, but the time one will be replaced with each index in `indices`
    # in the `for i in indices` loop below.
    idx = Any[Colon() for ind in axes(ses)]

    for i in indices
        idx[axis] = i
        namei = format("{}{:08d}", name, i)
        sesi = @view ses[idx...]
        vtkfile = to_vtp(namei, sesi)

        # Assume that all the times in sesi are the same. They should be if the
        # user specified the correct axis. I suppose I could check that and
        # issue a warning.
        t = first(sesi).τ
        pvd[t] = vtkfile
    end

    return pvd
end
