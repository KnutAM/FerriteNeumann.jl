module FerriteNeumann
using Ferrite

include("FerritePRs.jl")    # Ferrite PR's having required functionality
include("FerritePR495.jl")  # FaceIterator
include("defaultvalues.jl")
include("Neumann.jl")       
include("BodyLoad.jl")

export NeumannHandler, Neumann, BodyLoad

"""
    NeumannHandler(dh::AbstractDofHandler)    

The handler of all the Neumann boundary conditions in `dh`
that can be used to apply the neumann contribution to the 
external "force"-vector, `fext`:

```julia
fext = zeros(ndofs(dh))
nh=NeumannHandler(dh)
add!(nh, Neumann(...))  # Add boundary conditions
for t in timesteps
    fill!(fext, 0)
    apply!(fext, nh, t) # Add contributions to `fext`
    ...
end
```
"""
struct NeumannHandler{DH<:AbstractDofHandler}
    nbcs::Vector{NeumannData}
    bodyloads::Vector{BodyLoadData}
    dh::DH
end
NeumannHandler(dh::AbstractDofHandler) = NeumannHandler(NeumannData[], BodyLoadData[], dh)

Ferrite.add!(nh::NeumannHandler, nbc::Neumann) = add_neumann!(nh.nbcs, nbc, nh.dh)
Ferrite.add!(nh::NeumannHandler, bl::BodyLoad) = add_bodyload!(nh.bodyloads, bl, nh.dh)

# Application of boundary conditions
function Ferrite.apply!(f::Vector, nh::NeumannHandler, time)
    foreach(nbc->apply_neumann!(f,nbc,nh.dh,time), nh.nbcs)
    foreach(bl->apply_bodyload!(f,bl,nh.dh,time), nh.bodyloads)
end

include("show.jl")

end
