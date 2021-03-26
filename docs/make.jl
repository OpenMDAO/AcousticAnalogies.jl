using Documenter, AcousticAnalogies

IN_CI = get(ENV, "CI", nothing)=="true"

makedocs(sitename="AcousticAnalogies.jl", modules=[AcousticAnalogies], doctest=false,
         format=Documenter.HTML(prettyurls=IN_CI))

if IN_CI
    deploydocs(repo="github.com/dingraha/AcousticAnalogies.jl.git", devbranch="main")
end
