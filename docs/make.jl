using Documenter, FinEtools, FinEtoolsDeforLinear, FinEtoolsFlexBeams

makedocs(
	modules = [FinEtoolsFlexBeams],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsFlexBeams.jl",
	pages = Any[
			"Home" => "index.md",
			"How to guide" => "guide/guide.md",
			"Reference" => "man/reference.md"	
		],
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsFlexBeams.jl.git",
)
