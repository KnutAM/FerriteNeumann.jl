var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FerriteNeumann","category":"page"},{"location":"#FerriteNeumann","page":"Home","title":"FerriteNeumann","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FerriteNeumann.jl:  Simple Neumann boundary conditions for Ferrite.jl","category":"page"},{"location":"","page":"Home","title":"Home","text":"NeumannHandler\nNeumann","category":"page"},{"location":"#FerriteNeumann.NeumannHandler","page":"Home","title":"FerriteNeumann.NeumannHandler","text":"NeumannHandler(dh::AbstractDofHandler)\n\nThe handler of all the Neumann boundary conditions in dh.\n\nWith nh=NeumannHandler(dh), add boundary conditions  with add!(nh, Neumann(...)). \n\nApply the boundary conditions for a given time to the  external \"force\"-vector fext as apply!(fext, nh, time)\n\n\n\n\n\n","category":"type"},{"location":"#FerriteNeumann.Neumann","page":"Home","title":"FerriteNeumann.Neumann","text":"Neumann(field_name::Symbol, fv::FaceValues, faceset::Set{FaceIndex}, f, cellset=nothing)\n\nDefine a Neumann contribution with the weak forms according to \n\nint_Gamma f  delta u  mathrmdGamma quad textScalar fields \n\nint_Gamma boldsymbolf cdot boldsymboldelta u  mathrmdGamma\nquad textVector fields \n\nwhere Gamma is the boundary where the contribution is active.  f, or boldsymbolf, is the prescribed Neumann value,  defined by a function with signatures\n\nf(x::Vec, time, n::Vec) -> Number (Scalar field)\n\nf(x::Vec, time, n::Vec) -> Vec{dim} (Vector field)\n\nwhere x is the spatial position of the current quadrature point, time is the  current time, and n is the face normal vector. \n\nfield_name describes the field on which the boundary condition should abstract\nfv describes the interpolation and integration rules. \nfaceset describes which faces the BC is applied to\nif cellset!=nothing, only the cells in faceset that are also in cellset are included  (important when using mixed grids with different cells)\n\n\n\n\n\n","category":"type"}]
}
