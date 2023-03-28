# FerriteNeumann

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteNeumann.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteNeumann.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteNeumann.jl)

Simplified application of external loads, such as Neumann boundary conditions and body loads for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/).
Supports both `DofHandler` and `MixedDofHandler`, 
and can automatically infer suitable `FaceValues`/`CellValues` 
for a given the quadrature order. 
It is also possible to specify `FaceValues`/`CellValues` manually. 

### Neumann BCs

- Scalar field: $\int_{\Gamma} f \ \delta u \ \mathrm{d}\Gamma$
- Vector field: $\int_{\Gamma} \boldsymbol{f} \cdot \boldsymbol{\delta u} \ \mathrm{d}\Gamma$

### Body load

- Scalar field: $\int_{\Omega} f \ \delta u \ \mathrm{d}\Omega$
- Vector field: $\int_{\Omega} \boldsymbol{f} \cdot \boldsymbol{\delta u} \ \mathrm{d}\Omega$


## Setting up Neumann BC
Let's consider the case of a coupled problem with one 
scalar, `:c`, and one vector, `:u`, field. 

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

## Setting up body loads
Let's again consider the case of a coupled problem with one 
scalar, `:c`, and one vector, `:u`, field.

*While only body loads are considered in this example,* 
*the intended usage is to combine NeumannBCs and body loads* 
*in the same `NeumannHandler`.*

We start by defining how the body loads should vary with 
position (within the given cell set) and time.
Let's say we want to apply centrifugal body/volume force with a ramping rotational speed $\omega(t)=kt$ around the origin, i.e. $\boldsymbol{f}=\rho\omega^2\boldsymbol{r}$ on 
the displacement field, `:u`. 
On the concentration field, `:c`, we add a constant source term 
for all cells in the cellset `"supply"`. 
```julia
k = 1               # 1/s²
x0 = zero(Vec{3})   
f_u(x::Vec, time) = (r=x-x0; ω=k*t; (ρ*ω^2)*r) # ::Vec
f_c(x::Vec, time) = 1.0    # kg/(m³s) (::Number)
```

With the time and spatial variations defined, we can setup the actual 
loading types,
```julia
nh = NeumannHandler(dh::DofHandler)      # Container for Neumann BCs and body loads
add!(nh, BodyLoad(:c, 2, getcellset(grid, "supply"), f_c))   # Use 2nd order default quadrature integration
add!(nh, BodyLoad(:u, 2, f_u))   # Use 2nd order default quadrature integration and the entire domain (since no cellset is given)
```
In this case, the choice of `CellScalarValues` or `CellVectorValues` 
is made based on the return type of `f_c` and `f_u` in each case. 
The reference shape, function interpolation and 
default geometry interpolation is taken from the `dh`. 
The same code can be used for `DofHandler` and `MixedDofHandler`. 
In the same way as for `Neumann`, `CellValues` or `QuadratureRule`can 
be specified manually, with the same requirements regarding mixed grids. 

### Application
During the time stepping, the current forces can be added to the force vector with 

```julia
apply!(f::Vector, nh, time)
```
in each time step. Before this, `f` must be zeroed; either via `fill!(f, 0)` or with `Ferrite.start_assemble` if `f` is included in the assembler.
