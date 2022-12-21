"""
    Neumann(field_name::Symbol, facevalues::FaceValues, faceset::Set{FaceIndex}, f, cellset=nothing)

Define a Neumann contribution with the weak forms according to 
```math
\\int_{\\Gamma} f \\ \\delta u \\ \\mathrm{d}\\Gamma, \\quad \\text{Scalar fields} \\\\

\\int_{\\Gamma} \\boldsymbol{f} \\cdot \\boldsymbol{\\delta u} \\ \\mathrm{d}\\Gamma,
\\quad \\text{Vector fields} 
```
where ``\\Gamma`` is the boundary where the contribution is active. 
``f``, or ``\\boldsymbol{f}``, is the prescribed Neumann value, 
defined by a function with signatures

`f(x::Vec, time, n::Vec) -> Number` (Scalar field)

`f(x::Vec, time, n::Vec) -> Vec{dim}` (Vector field)

where `x` is the spatial position of the current quadrature point, `time` is the 
current time, and `n` is the face normal vector. 

* `fieldname` describes the field on which the boundary condition should abstract
* `facevalues` describes the interpolation and integration rules. 
* `faceset` describes which faces the BC is applied to
* if `cellset!=nothing`, only the cells in `faceset` that are also in `cellset` are included 
  (important when using mixed grids with different cells)
"""
struct Neumann{FV,FUN}
    fieldname::Symbol
    facevalues::FV
    faceset::Set{FaceIndex}
    f::FUN # f(x::Vec, time, n::Vec)->{FV::FaceScalarValues ? Number : Vec}
end

function Neumann(fieldname::Symbol, facevalues::FaceValues, faceset::Set{FaceIndex}, f, cellset)
    _faceset = intersect_faces_and_cells(faceset, cellset)
    return Neumann(fieldname, facevalues, _faceset, f)
end

intersect_faces_and_cells(faceset::Set{FaceIndex}, ::Nothing) = faceset
function intersect_faces_and_cells(faceset::Set{FaceIndex}, cellset)
    return Set(face for face in faceset if first(face) in cellset)
end


struct NeumannHandler{DH<:AbstractDofHandler}
    nbcs::Vector
    dh::DH
end

"""
    NeumannHandler(dh::AbstractDofHandler)    

The handler of all the Neumann boundary conditions in `dh`.

With `nh=NeumannHandler(dh)`, add boundary conditions 
with `add!(nh, Neumann(...))`. 

Apply the boundary conditions for a given `time` to the 
external "force"-vector `fext` as `apply!(fext, nh, time)`
"""
NeumannHandler(dh::AbstractDofHandler) = NeumannHandler([], dh)

function Ferrite.add!(nh::NeumannHandler, nbc)
    # Should make some checks here...
    push!(nh.nbcs, nbc)
end

function Ferrite.apply!(f::Vector, nh::NeumannHandler, time)
    foreach(nbc->apply!(f,nbc,nh.dh,time), nh.nbcs)
end

function Ferrite.apply!(f::Vector{T}, nbc::Neumann, dh::DofHandler, time) where T
    dofs = collect(dof_range(dh, nbc.fieldname))
    fe = zeros(T, length(dofs))
    for face in FaceIterator(dh, nbc.faceset)
        calculate_neumann_contribution!(fe, face, nbc.facevalues, time, nbc.f)
        assemble!(f, view(celldofs(face), dofs), fe)
    end
end

function calculate_neumann_contribution!(fe::Vector, face::FaceCache, fv::FaceValues, time, f)
    fill!(fe, 0)
    reinit!(fv, face)
    for q_point in 1:getnquadpoints(fv)
        dΓ = getdetJdV(fv, q_point)
        x = spatial_coordinate(fv, q_point, getcoordinates(face))
        n = getnormal(fv, q_point)
        b = f(x, time, n)
        for i in 1:getnbasefunctions(fv)
            δu = shape_value(fv, q_point, i)
            fe[i] += (δu ⋅ b) * dΓ
        end
    end
end