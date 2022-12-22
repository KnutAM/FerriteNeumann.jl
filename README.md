# FerriteNeumann

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteNeumann.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteNeumann.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteNeumann.jl)

Simplified application of Neumann boundary conditions for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/), 
where the Neumann BCs are of the type

- Scalar field: $\int_{\Gamma} f \ \delta u \ \mathrm{d}\Gamma$
- Vector field: $\int_{\Gamma} \boldsymbol{f} \cdot \boldsymbol{\delta u} \ \mathrm{d}\Gamma$

Works conveniently for both `DofHandler` and `MixedDofHandler` and can automatically create the 
suitable `FaceValues` when given the quadrature order. It is also possible to specify `FaceValues`
manually. 

## Example
Let's consider the case of a coupled problem with one 
scalar, `:c`, and one vector, `:u`, field. 

### Definition
We start by defining how the Neumann boundary conditions should vary with 
position (within the given face set), time, and normal vector. 
Let's say we want to apply a pressure as a function of height (`x[3]`) on 
the displacement field, `:u`, and a ramping diffusion flux for the 
concentration field, `:c`. 
```julia
ρg = 1e4    # N/m³
f_u(x::Vec, time, n::Vec) = -x[3]*n # ::Vec, normal to face (in this case)
f_c(x::Vec, time, n::Vec) = time    # ::Number
```

With the time and spatial variation defined, we can setup the actual 
boundary value types,
```julia
nh = NeumannHandler(dh::DofHandler)      # Container for Neumann BCs
faceset = getfaceset(grid, "right")      # 
add!(nh, Neumann(:c, 2, faceset, f_c))   # Use 2nd order default quadrature integration 
add!(nh, Neumann(:u, 2, faceset, f_u))   # Use 2nd order default quadrature integration 
```
In this case, the choice of `FaceScalarValues` or `FaceVectorValues` is made based 
on the return type of `f_c` and `f_u` in each case. The reference shape, function 
interpolation and default geometry interpolation is taken from the `dh`. 
The same code can be used for `DofHandler` and `MixedDofHandler`.

<details>
<summary>Specifying `FaceValues` or `QuadratureRule` manually</summary>
It is also possible to specify either the specific `QuadratureRule` manually, 
e.g. 

```julia
qr = QuadratureRule{dim-1, RefCube}(2)    # 
add!(nh, Neumann(:c, qr, faceset, f_c))   # 
add!(nh, Neumann(:u, qr, faceset, f_u))   # 
```

or even the complete `FaceValues`
```julia
ip = Lagrange{dim, RefCube, 1}()            #
fv_c = FaceScalarValues(qr, ip)             # :c requires scalar values
fv_u = FaceVectorValues(qr, ip)             # :u requires vector values
add!(nh, Neumann(:c, fv_c, faceset, f_c))   # Use 2nd order default quadrature integration 
add!(nh, Neumann(:u, fv_u, faceset, f_u))   # Use 2nd order default quadrature integration 
```

Note that in these cases, care must be taken to only have faces in `faceset` that 
are compatible with the given `FaceValues` or `QuadratureRule`. This may be an issue when 
using the `MixedDofHandler`. 

</details>  
&nbsp;  

### Application
During the time stepping, the current forces can be added to the force vector with 

```julia
apply!(f::Vector, nh, time)
```
noting that the full force is added (not only the increment), so this statement 
is usually prepended by `fill!(f, 0)`. However, this action can also be done before 
the element assembly if it includes `f`, nothing that it must be done after 
`start_assemble` as that zeros `f`.  
