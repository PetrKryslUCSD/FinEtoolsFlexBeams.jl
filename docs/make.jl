using Documenter, FinEtools, FinEtoolsLinearDefor, FinEtoolsFlexBeams

makedocs(
	modules = [FinEtoolsFlexBeams],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsFlexBeams.jl",
	pages = Any[
	"Home" => "index.md",
	"Tutorials" => "tutorials/tutorials.md",
	"Types and Functions" => Any[
		"man/types.md",
		"man/functions.md"]
		]
	"Guide" => "guide/guide.md",
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsFlexBeams.jl.git",
)
