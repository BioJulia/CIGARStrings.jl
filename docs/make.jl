using Documenter, CIGARStrings

meta = quote
    using CIGARStrings
    cigar = CIGAR("9H1S15M1D3X3=17M12I8S")
end

DocMeta.setdocmeta!(CIGARStrings, :DocTestSetup, meta; recursive=true)

makedocs(;
    sitename="CIGARStrings.jl",
    modules=[CIGARStrings],
    pages=["Home" => "index.md"],
    authors="Jakob Nybo Nissen",
    checkdocs=:public,
    remotes=nothing,
)

deploydocs(;
    repo="github.com/BioJulia/CIGARStrings.jl.git",
    push_preview=true,
    deps=nothing,
    make=nothing,
)
