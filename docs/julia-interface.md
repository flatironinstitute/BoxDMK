# Julia Interface

This package exposes:

- low-level wrappers around the Fortran/C ABI
- a small high-level workflow API for common solve paths

## Low-level API

Low-level wrappers are:

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

### Callback contract for tree builders

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

### Required buffer shapes

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

## High-level API

High-level types/functions are implemented in `julia/src/highlevel.jl`:

- `BDMKProblem`
- `BDMKOptions`
- `BDMKTree`
- `BDMKResult`
- `build_tree`
- `solve`
- `evaluate_targets`
- `run`

### `BDMKProblem` parameters

`BDMKProblem` is created with keyword arguments:

- `density::Function`: callback called as `density(x, problem)` where `x` is a `Vector{Float64}` of length `ndim`. Return either one number (broadcast to all components) or a vector of length `nd`.
- `nd::Integer` (default `1`): number of components in field values.
- `ndim::Integer` (default `3`): spatial dimension (`1`, `2`, or `3`).
- `ikernel::Integer` (default `1`): kernel selector passed to `bdmk!`.
- `beta::Real` (default `1.0`): kernel parameter passed to `bdmk!`.
- `boxlen::Real` (default `1.18`): computational box side length.
- `dpars::Vector{Float64}`: real parameter buffer used by callback/kernel.
- `zpars::Vector{ComplexF64}`: complex parameter buffer used by callback/kernel.
- `ipars::Vector{Cint}`: integer parameter buffer used by callback/kernel.

The constructor normalizes key entries in `ipars`:

- `ipars[1] = ndim`
- `ipars[2] = ikernel`

### `BDMKOptions` parameters

- `eps::Float64` (default `1e-6`): requested solve tolerance passed to `bdmk!`.
- `norder::Int` (default `16`): polynomial order per axis; points per box are `npbox = norder^ndim`.
- `ipoly::Int` (default `0`): polynomial type flag for tree construction.
- `iptype::Int` (default `2`): tree node type flag.
- `eta::Float64` (default `0.0`): tree builder parameter passed through to low-level calls.
- `epstree_factor::Float64` (default `500.0`): tree tolerance scaling (`eps_tree = eps * epstree_factor`).
- `ifnewtree::Int` (default `0`): tree-construction mode flag.
- `iperiod::Int` (default `0`): periodicity flag.
- `zk::ComplexF64` (default `30.0 + 0.0im`): complex kernel parameter for tree calls.

### Workflow functions

- `tree = build_tree(problem, opts=BDMKOptions())`: performs the low-level two-stage tree build (`vol_tree_mem!`, then `vol_tree_build!`) and returns `BDMKTree`.
- `result = solve(problem, tree; compute=:potential, targets=nothing, eps=1e-6)`: evaluates on tree nodes, and optionally on target points.
- `result = evaluate_targets(problem, tree, targets; compute=:potential, eps=1e-6)`: convenience wrapper around `solve(...; targets=...)`.
- `(tree, result) = run(problem; targets=nothing, compute=:potential, opts=BDMKOptions())`: one-call build + solve.

`compute` accepts:

- `:potential` (`ifpgh=1`)
- `:gradient` (`ifpgh=2`)
- `:hessian` (`ifpgh=3`)

`targets` must have shape `(ndim, ntarg)`.

### Return types

`BDMKTree` contains all tree/buffer data needed by `bdmk!` (`itree`, `iptr`, `centers`, `boxsize`, `fvals`, metadata).

`BDMKResult` fields:

- `pot`: node values, shape `(nd, npbox, nboxes)`
- `grad`: node gradients or `nothing`
- `hess`: node Hessians or `nothing`
- `pote`: target values or `nothing`
- `grade`: target gradients or `nothing`
- `hesse`: target Hessians or `nothing`
- `tottimeinfo`: timing array from low-level solver
- `ifpgh`, `ifpghtarg`: internal output level flags

### Quick start

```julia
include("julia/src/BoxDMK.jl")
using .BoxDMK

function density(x, problem)
    r2 = sum(abs2, x)
    return exp(-10r2)
end

prob = BDMKProblem(density=density, nd=1, ndim=3, ikernel=1)
opts = BDMKOptions(eps=1e-6, norder=8)

tree, res = BoxDMK.run(prob; compute=:potential, opts=opts)

targets = rand(3, 20) .- 0.5
tres = evaluate_targets(prob, tree, targets; compute=:potential, eps=opts.eps)
```

### Thread-safety caveat

The high-level callback trampoline uses a single global callback context. Concurrent tree builds from multiple Julia threads/tasks are not safe with the current implementation.
