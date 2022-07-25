"""
    Neumann(field_name::Symbol, fv::FaceValues, faceset::Set{FaceIndex}, f, cellset=nothing)

    Define a Neumann contribution according to 
    ```math
    \\int_{\\Gamma} b \\ \\delta u \\ \\mathrm{d}\\Gamma

    \\int_{\\Gamma} \\boldsymbol{b} \\cdot \\boldsymbol{\\delta u} \\ \\mathrm{d}\\Gamma
    ```
    for scalar or vector fields, respectively. Here ``b`` or 
    ``\\boldsymbol{b}`` is the prescribed Neumann value, 
    defined by the function `f` with signature 
    ``` 
    f(x::Vec, time, n::Vec)
    ```
    where `x` is the spatial position of the current quadrature point, `time` is the 
    current time, and `n` is the face normal vector. `field_name` describes the 
    field on which the boundary condition should act, and `fv` describes the 
    interpolation and integration rules. `faceset` describes which faces the BC is 
    applied to, and if `cellset!=nothing`, only the cells in `faceset` that are also 
    in `cellset` are included (important when using mixed grids with different cells)
"""
struct Neumann{FV,FUN}
    field_name::Symbol
    face_values::FV
    faces::Vector{Int}
    cells::Vector{Int}
    f::FUN # f(fv::FaceValues, q_point::Int, i_shape::Int, x::Vec, t::Number)->Number
end

function Neumann(field_name::Symbol, fv::FaceValues, faceset::Set{FaceIndex}, f, cellset=nothing)
    cells, faces = _get_cells_and_faces(faceset, cellset)
    return Neumann(field_name, fv, faces, cells, f)
end

struct NeumannHandler{DH<:AbstractDofHandler}
    nbcs::Vector
    dh::DH
end

NeumannHandler(dh::AbstractDofHandler) = NeumannHandler([], dh)

function Ferrite.add!(nh::NeumannHandler, nbc)
    # Should make some checks here...
    push!(nh.nbcs, nbc)
end

function Ferrite.apply!(f::Vector, nh::NeumannHandler, time)
    foreach(nbc->apply!(f,nbc,nh.dh,time), nh.nbcs)
end

function Ferrite.apply!(f::Vector{T}, nbc::Neumann, dh::DofHandler, time) where T
    dofs = collect(dof_range(dh, nbc.field_name))
    fe = zeros(T, length(dofs))
    for face in FaceIterator(dh, nbc.faces, nbc.cells)
        calculate_neumann_contribution!(fe, face, nbc.face_values, time, nbc.f)
        assemble!(f, view(celldofs(face), dofs), fe)
    end
end

function calculate_neumann_contribution!(fe::Vector, face::FaceIterator, fv::FaceValues, time, f)
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