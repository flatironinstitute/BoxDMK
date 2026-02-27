# Julia Low-Level Interface

This interface exposes three low-level routines through `ccall`:

- `vol_tree_mem!`
- `vol_tree_build!`
- `bdmk!`

The wrappers are implemented in `julia/src/BoxDMK.jl` and call C-ABI Fortran shims exported from `libboxdmk.so`.

## Build shared library

```bash
module purge
module load gcc cmake lib/openblas
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Expected artifact: `build/libboxdmk.so`.

## Run Julia tests

```bash
julia --project=julia julia/test/runtests.jl
```

## Load module

```julia
include("julia/src/BoxDMK.jl")
using .BoxDMK
```

By default, wrappers load `build/libboxdmk.so` relative to repo root.
You can override this:

```julia
BoxDMK.set_libboxdmk_path!("/absolute/path/to/libboxdmk.so")
```

## Callback contract for tree builders

`vol_tree_mem!` and `vol_tree_build!` require a callback pointer (`funptr`) created by `@cfunction`.

Expected callback ABI (C style):

```julia
# Conceptual signature
# nd is passed by value, other arrays are pointers.
Cvoid callback(Cint nd,
               Ptr{Cdouble} xyz,
               Ptr{Cdouble} dpars,
               Ptr{ComplexF64} zpars,
               Ptr{Cint} ipars,
               Ptr{Cdouble} f)
```

The callback must write `nd` output values into `f`.

## Required buffer shapes

Caller is responsible for allocating all buffers.

- `rintl`: length >= 201 for `vol_tree_mem!`
- `iptr`: length >= 8
- `itree`: length `ltree`
- `centers`: length `ndim*nboxes`
- `boxsize`: length `nlevels+1`
- `fvals`: length `nd*(norder^ndim)*nboxes`

For `bdmk!`, typical flattened lengths are:

- `pot`: `nd*npbox*nboxes`
- `grad`: `nd*ndim*npbox*nboxes` (when `ifpgh >= 2`)
- `hess`: `nd*(ndim*(ndim+1)/2)*npbox*nboxes` (when `ifpgh >= 3`)
- `targs`: `ndim*ntarg`
- `pote`: `nd*ntarg`
- `grade`: `nd*ndim*ntarg`
- `hesse`: `nd*(ndim*(ndim+1)/2)*ntarg`

All arrays are Fortran-contiguous logically; Julia vectors should be flattened in column-major order.
