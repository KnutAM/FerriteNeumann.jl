import Ferrite: ScalarWrapper, AbstractDofHandler   # Used, but not exported by Ferrite.jl
import Ferrite: reinit!, celldofs!, onboundary      # Overloaded
#= # Only if merged with Ferrite.jl, currently no AbstractGridIterator
abstract type AbstractGridIterator end
Base.IteratorSize(::Type{T})   where {T<:AbstractGridIterator} = Base.HasLength() # this is default in Base
Base.IteratorEltype(::Type{T}) where {T<:AbstractGridIterator} = Base.HasEltype() # this is default in Base
Base.eltype(::Type{T})         where {T<:AbstractGridIterator} = T
=#


struct FaceIterator{CI<:CellIterator}# <: AbstractGridIterator
    faces::Vector{Int}
    current_faceid::ScalarWrapper{Int}
    ci::CI
end

function _get_cells_and_faces(faceset, ::Nothing)
    tuple((collect([faceindex[j] for faceindex in faceset]) for j in 1:2)...)
end

function _get_cells_and_faces(faceset, cellset)
    tuple((collect([faceindex[j] for faceindex in faceset if faceindex[1] in cellset]) for j in 1:2)...)
end

function FaceIterator(dh::AbstractDofHandler, faceset::Set{FaceIndex}, cellset=nothing, args...)
    cells, faces = _get_cells_and_faces(faceset, cellset)
    return FaceIterator(faces, ScalarWrapper(0), CellIterator(dh, cells, args...))
end

function FaceIterator(dh::AbstractDofHandler, faces::Vector{Int}, cells::Vector{Int}, args...)
    if length(faces)!=length(cells)
        msg = "faces and cells have different lengths: $(length(faces)) vs $(length(cells))"
        throw(DimensionMismatch(msg))
    end
    return FaceIterator(faces, ScalarWrapper(0), CellIterator(dh, cells, args...))
end

# Can be removed an replaced by AbstractGridIterator functions if internal in Ferrite.jl
Base.IteratorSize(::Type{T})   where {T<:FaceIterator} = Base.HasLength() # this is default in Base
Base.IteratorEltype(::Type{T}) where {T<:FaceIterator} = Base.HasEltype() # this is default in Base
Base.eltype(::Type{T})         where {T<:FaceIterator} = T

@inline Base.length(fi::FaceIterator)  = length(fi.ci)
function Base.iterate(fi::FaceIterator, state = 1)
    if state > length(fi)
        return nothing
    else
        return (reinit!(fi, state), state+1)
    end
end

# Use functions from CellIterator, except nfaces
# New functions: faceid, faceindex
function reinit!(fi::FaceIterator, i::Int)
    reinit!(fi.ci, i)
    fi.current_faceid[] = fi.faces[i]
    return fi
end

for op = (:getnodes, :getcoordinates, :cellid, :celldofs)
    eval(quote
        function Ferrite.$op(fi::FaceIterator, args...; kwargs...)
            return Ferrite.$op(fi.ci, args...; kwargs...)
        end
    end)
end
@inline faceid(fi::FaceIterator) = fi.current_faceid[]
@inline celldofs!(v::Vector, fi::FaceIterator) = celldofs!(v, fi.ci)
@inline onboundary(fi::FaceIterator) = onboundary(fi.ci, faceid(fi))
@inline faceindex(fi::FaceIterator) = FaceIndex(cellid(fi), faceid(fi))
@inline function reinit!(fv::FaceValues, fi::FaceIterator)
    reinit!(fv, fi.ci, faceid(fi))
end