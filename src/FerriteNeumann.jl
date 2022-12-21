module FerriteNeumann
using Ferrite

include("PR495.jl")     # FaceIterator
include("neumann.jl")

export NeumannHandler, Neumann

end
