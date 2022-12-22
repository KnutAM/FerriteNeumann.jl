"""
    Neumann(field_name::Symbol, fv_info::Union{FaceValues,QuadratureRule,Int}, faceset::Set{FaceIndex}, f)

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
current time, and `n` is the face normal vector. The remaining input arguments are

* `fieldname` describes the field on which the boundary condition should abstract
* `fv_info` gives required input to determine the facevalues. The following input types are accepted:
  - `FaceValues` matching the interpolation for `fieldname` for the faces in `faceset`
  - `QuadratureRule` matching the interpolation for `fieldname` for faces in `faceset`
  - `Int` giving the integration order to use. `FaceValues` are deduced from the interpolation 
    of `fieldname` and the output of `f`. 
* `faceset` describes which faces the BC is applied to
"""
struct Neumann{FVI,FUN}
    fieldname::Symbol
    fv_info::FVI
    faceset::Set{FaceIndex}
    f::FUN # f(x::Vec, time, n::Vec)->{FV::FaceScalarValues ? Number : Vec}
end

struct NeumannData{FV,FUN}
    fieldname::Symbol   # Only for information 
    dofrange::UnitRange{Int}
    facevalues::FV 
    faceset::Set{FaceIndex}
    f::FUN
end

function NeumannData(dh::DofHandler, spec::Neumann)
    dofrange = dof_range(dh, spec.fieldname)
    Ferrite.get_func_interpolations(dh, spec.fieldname)
    ip = Ferrite.getfieldinterpolation(dh, Ferrite.find_field(dh, spec.fieldname))
    fv = get_facevalues(spec.fv_info, ip, spec.f)
    return NeumannData(spec.fieldname, dofrange, fv, spec.faceset, spec.f)
end

# If a facevalue has already been given, use this value 
get_facevalues(fv::FaceValues, args...) = fv

# Create default quadrule of given order
function get_facevalues(order::Int, ip::Interpolation{dim,RefShape}, f) where {dim, RefShape}
    return get_facevalues(QuadratureRule{dim-1,RefShape}(order), ip, f)
end

# Use the given function to determine if the output should be a scalar or vector. 
get_facevalues(qr::QuadratureRule, ip::Interpolation{dim}, f) where dim = get_facevalues(qr, ip, f(zero(Vec{dim}), 0.0, zero(Vec{dim})))
function get_facevalues(qr::QuadratureRule{<:Any,RefShape}, ip::Interpolation{<:Any,RefShape}, ::Vec) where RefShape
    return FaceVectorValues(qr, ip)
end
function get_facevalues(qr::QuadratureRule{<:Any,RefShape}, ip::Interpolation{<:Any,RefShape}, ::Number) where RefShape
    return FaceScalarValues(qr, ip)
end
function get_facevalues(qr::QuadratureRule, ip::Interpolation, fval::Union{Vec,Number})
    throw(ArgumentError("qr, $(typeof(qr)), and ip, $(typeof(ip)), doesn't seem compatible. (info: fval=$fval)"))
end

struct NeumannHandler{DH<:AbstractDofHandler}
    nbcs::Vector
    dh::DH
end

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
NeumannHandler(dh::AbstractDofHandler) = NeumannHandler([], dh)

function Ferrite.add!(nh::NeumannHandler{<:DofHandler}, nbc::Neumann)
    push!(nh.nbcs, NeumannData(nh.dh, nbc))
end

function Ferrite.add!(nh::NeumannHandler{<:MixedDofHandler}, nbc::Neumann)
    for fh in nh.dh.fieldhandlers
        push!(nh.nbcs, NeumannData(nh.dh, fh, nbc))
    end
end

function Ferrite.apply!(f::Vector, nh::NeumannHandler, time)
    foreach(nbc->apply!(f,nbc,nh.dh,time), nh.nbcs)
end

function Ferrite.apply!(f::Vector{T}, nbc::NeumannData, dh::DofHandler, time) where T
    dofs = collect(nbc.dofrange)
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
