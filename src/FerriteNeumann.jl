module FerriteNeumann
using Ferrite

include("FerritePRs.jl")    # Old PR's having required functionality
include("PR495.jl")         # FaceIterator
include("neumann.jl")       

export NeumannHandler, Neumann

end
