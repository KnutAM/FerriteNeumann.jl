# PR648: hasfieldname
@static if hasproperty(Ferrite, :hasfieldname)
    import Ferrite: hasfieldname
else
    hasfieldname(fh::FieldHandler, field_name::Symbol) = !isnothing(findfirst(field->field.name==field_name, fh.fields))
    #hasfieldname(dh::Ferrite.AbstractDofHandler, field_name::Symbol) = field_name ∈ Ferrite.getfieldnames(dh)
end