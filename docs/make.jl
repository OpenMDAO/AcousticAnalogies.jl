module AADocs
using Documenter, AcousticAnalogies

function doit()
    IN_CI = get(ENV, "CI", nothing)=="true"

    makedocs(sitename="AcousticAnalogies.jl", modules=[AcousticAnalogies], doctest=false,
             root=@__DIR__,
             format=Documenter.HTML(prettyurls=IN_CI),
             pages=["Introduction"=>"index.md",
                    "Guided Example"=>"guided_example.md",
                    "CCBlade.jl Example"=>"ccblade_example.md",
                    "WriteVTK.jl Support"=>"writevtk_support.md",
                    "OpenFAST Example"=>"openfast_example.md",
                    "API Reference"=>"api.md",
                    "Software Quality Assurance"=>"sqa.md",
                    "BPM Airfoil Self-Noise Tests"=>"bpm_tests1.md",
                    "BPM Airfoil Self-Noise Tests, Cont."=>"bpm_tests2.md",
                    "BPM Airfoil Self-Noise Tests, Cont."=>"bpm_tests3.md",
                    "Ideally Twisted Rotor Tests"=>"itr_tests1.md",
                    "Ideally Twisted Rotor Tests, Cont."=>"itr_tests2.md",
                   ])

    if IN_CI
        deploydocs(repo="github.com/OpenMDAO/AcousticAnalogies.jl.git", devbranch="main")
    end

end

if !isinteractive()
    doit()
end

end # module
