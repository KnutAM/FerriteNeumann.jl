# FerriteNeumann

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://KnutAM.github.io/FerriteNeumann.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://KnutAM.github.io/FerriteNeumann.jl/dev)
[![Build Status](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/KnutAM/FerriteNeumann.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/KnutAM/FerriteNeumann.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/KnutAM/FerriteNeumann.jl)

Simplified application of Neumann boundary conditions for [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl/), 
where the Neumann BCs are of the type

- Scalar field: $\int_{\Gamma} b \ \delta u \ \mathrm{d}\Gamma$
- Vector field: $\int_{\Gamma} \boldsymbol{b} \cdot \boldsymbol{\delta u} \ \mathrm{d}\Gamma$

## Example
Let's consider the case of a coupled problem with one 
scalar, `:c`, and one vector, `:u`, field. 

We start by defining how the Neumann boundary conditions should vary with 
position (within the given face set), time, and normal vector. 
Let's say we want to apply a pressure as a function of height (`x[3]`) on 
the displacement field, `:u`, and a ramping diffusion flux for the 
concentration field, `:c`. 
```julia
function f_u(x::Vec, time, n::Vec)
    ρ = 1000; g = 9.81
    return -ρ*g*x[3]*n              # Returns a Vec
end
f_c(x::Vec, time, n::Vec) = time    # Returns a Number
```

Note that for scalar fields, `fv::FaceScalarValues` and `f` should 
return a scalar. For vector fields, `fv::FaceVectorValues` and `f` 
should return a `Vec` with the same dimension as the field.

With the time and spatial variation defined, we can setup the actual 
boundary value types, starting with the container, defining the faceset, 
the facevalues, and finally adding each condition

```julia
nh = NeumannHandler(dh::AbstractDofHandler) # Container for Neumann BCs
faceset = getfaceset(grid, "right")         # 
qr = QuadratureRule{dim-1, RefCube}(2)      # 
ip = Lagrange{dim, RefCube, 1}()            #
fv_c = FaceScalarValues(qr, ip)             # :c requires scalar values
fv_u = FaceVectorValues(qr, ip)             # :u requires vector values
add!(nh, Neumann(:c, fv_c, faceset, f_c))   #
add!(nh, Neumann(:u, fv_u, faceset, f_u))   # 
```

During the time stepping, the current forces can be added to the force vector with 

```julia
apply!(f::Vector, nh, time)
```
noting that the full force vector is added (not only the increment), so this statement 
is usually prepended by `fill!(f, 0)`. However, this action can also be done before 
the element assembly if it includes `f`, nothing that start_assemble has this effect
if used on `f`. 