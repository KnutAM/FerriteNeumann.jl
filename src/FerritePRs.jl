# PR532: apply_analytical!
@static if !hasmethod(Ferrite.getfieldinterpolation, (FieldHandler, Int))
    Ferrite.getfieldinterpolation(fh::FieldHandler, field_idx::Int) = fh.fields[field_idx].interpolation
    Ferrite.getfielddim(fh::FieldHandler, field_idx::Int) = fh.fields[field_idx].dim
end