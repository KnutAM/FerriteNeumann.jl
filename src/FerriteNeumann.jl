module FerriteNeumann
using Ferrite

# Introduce FaceIterator, might be merged with Ferrite.jl
include("iterators.jl") 
export FaceIterator 
export faceid, faceindex

include("neumann.jl")   
export NeumannHandler, Neumann

end
