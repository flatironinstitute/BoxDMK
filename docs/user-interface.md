# User Interface Guide

This project is a Fortran library + test executable. It does not provide an interactive CLI.  
The practical user interface is:

- build/run commands (`make` or `cmake` + `ctest`)
- the public Fortran routine `subroutine bdmk(...)`
- the tree-construction routines used before calling `bdmk`

## 1) Runtime entrypoint

The executable `int2-bdmk` is built from `test/bdmk/testbdmk.f` and runs an end-to-end example:

- build adaptive tree data (`vol_tree_mem`, `vol_tree_build`)
- call `bdmk(...)`
- print timing, throughput, and relative error

Reference:
- `test/bdmk/testbdmk.f` calls `bdmk` near line 414.

## 2) Core API: `subroutine bdmk`

Defined in `src/bdmk/bdmk4.f`.

Signature (conceptual):

```fortran
subroutine bdmk(
  nd, ndim, eps, ikernel, beta, ipoly, norder, npbox,
  nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals,
  ifpgh, pot, grad, hess, ntarg, targs,
  ifpghtarg, pote, grade, hesse, tottimeinfo
)
```

### Required input controls

- `nd`: number of right-hand sides (densities).
- `ndim`: problem dimension (2 or 3 in current usage).
- `eps`: requested precision.
- `ikernel`: kernel selector.
  - `0`: Yukawa
  - `1`: Laplace
  - `2`: square-root Laplace
- `beta`: kernel parameter (Yukawa/power-kernel parameter).
- `ipoly`: tensor basis selector.
  - `0`: Legendre
  - `1`: Chebyshev
- `norder`: polynomial order per dimension.
- `npbox`: points per box (`norder**ndim`).

### Tree and source data inputs

- `nboxes, nlevels, ltree, itree, iptr, centers, boxsize`: adaptive tree structure and geometry.
- `fvals(nd, npbox, nboxes)`: source values on tensor grids in leaf boxes.

The test driver builds these arrays with:

- `vol_tree_mem(...)`
- `vol_tree_build(...)`

### Output selection flags

- `ifpgh`: what to compute on tree nodes.
  - `1`: potential only
  - `2`: potential + gradient
  - `3`: potential + gradient + Hessian
- `ifpghtarg`: what to compute at user targets.
  - `0`: no target evaluation
  - `1`: potential
  - `2`: potential + gradient
  - `3`: potential + gradient + Hessian

### Outputs

- On tree grid:
  - `pot(nd, npbox, nboxes)`
  - `grad(nd, ndim, npbox, nboxes)` when `ifpgh >= 2`
  - `hess(nd, ndim*(ndim+1)/2, npbox, nboxes)` when `ifpgh >= 3`
- At targets:
  - `pote(nd, ntarg)`
  - `grade(nd, ndim, ntarg)`
  - `hesse(nd, ndim*(ndim+1)/2, ntarg)`
- Timing:
  - `tottimeinfo(*)`

## 3) Memory/layout expectations

- Arrays are column-major Fortran layout.
- Allocate all output arrays before calling `bdmk`.
- Initialize output arrays (especially `pot/grad/hess`) before call, matching the test pattern.
- Hessian component ordering from `bdmk4.f` comments:
  - 2D: `xx, xy, yy`
  - 3D: `xx, yy, zz, xy, xz, yz`

## 4) Typical call flow

1. Set model/control parameters (`ndim`, `ikernel`, `eps`, `norder`, flags).
2. Build adaptive tree and source values via `vol_tree_mem` + `vol_tree_build`.
3. Allocate `pot/grad/hess` and optional target arrays.
4. Call `bdmk(...)`.
5. Postprocess: norms, error checks, throughput.

The reference implementation for this flow is `test/bdmk/testbdmk.f`.

## 5) Build-level user interface

Use one of:

- Makefile path (see `README.md`):
  - `make clean`
  - `make HOST=linux-gfortran FEND='-lopenblas' all`
- CMake path (see `README.md`):
  - load modules
  - configure/build/test with `cmake` and `ctest`

## 6) Known behavior from current code

- The default executable writes result logs like `l3dbox.txt`.
- Debug/progress prints are produced by `prini/prinf/prin2` routines.
- Current CMake/test build passes on OpenBLAS; one Fortran argument-type warning is present in the test code and does not fail the build.


## Tree Builder API

Tree construction is a two-step API in `src/common/tree_vol_coeffs.f`.

`vol_tree_mem` (memory sizing):

```fortran
subroutine vol_tree_mem(ndim,ipoly,iperiod,eps,zk,boxlen,norder,
1    iptype,eta,fun,nd,dpars,zpars,ipars,ifnewtree,nboxes,nlevels,
2    ltree,rintl)
```

`vol_tree_build` (actual build):

```fortran
subroutine vol_tree_build(ndim,ipoly,iperiod,eps,zk,boxlen,
1    norder,iptype,eta,fun,nd,dpars,zpars,ipars,rintl,nboxes,
2    nlevels,ltree,itree,iptr,centers,boxsize,fvals)
```

Parameter groups:
- Controls: `ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta`
- Density callback: `fun`
- Callback params: `nd, dpars, zpars, ipars`

Outputs:
- `vol_tree_mem`: `nboxes, nlevels, ltree, rintl`
- `vol_tree_build`: `itree, iptr, centers, boxsize, fvals`

### `vol_tree_mem` parameters (detailed)

Purpose: estimate tree size and memory requirements before allocation/build.

- `ndim` (`integer`, in): spatial dimension.
- `ipoly` (`integer`, in): polynomial family for node/grid generation.
  - `0`: Legendre
  - `1`: Chebyshev
- `iperiod` (`integer`, in): periodicity mode used by tree routines.
  - `0`: free-space (used in current tests)
  - `1`: periodic mode
- `eps` (`real*8`, in): requested tree/adaptivity accuracy.
- `zk` (`complex*16`, in): kernel parameter forwarded through tree setup code.
- `boxlen` (`real*8`, in): side length of the root box; domain is centered at origin.
- `norder` (`integer`, in): 1D polynomial order per box dimension.
- `iptype` (`integer`, in): error norm for refinement decisions.
  - `0`: Linf
  - `1`: L1
  - `2`: L2
- `eta` (`real*8`, in): refinement/error scaling factor.
- `fun` (`external`, in): user callback evaluating the source function values.
- `nd` (`integer`, in): number of scalar fields (RHS channels) returned by `fun`.
- `dpars(*)` (`real*8`, in): user real parameters passed to `fun`.
- `zpars(*)` (`complex*16`, in): user complex parameters passed to `fun`.
- `ipars(*)` (`integer`, in): user integer parameters passed to `fun`.
- `ifnewtree` (`integer`, in): memory policy hint.
  - `0`: assume unchanged tree
  - `1`: reserve extra memory for tree changes/refinement
- `nboxes` (`integer`, out): estimated total number of boxes.
- `nlevels` (`integer`, out): estimated number of levels.
- `ltree` (`integer`, out): required integer workspace length for packed tree arrays.
- `rintl(0:...)` (`real*8`, out): per-level norm estimates used for consistent error scaling; caller supplies enough storage (code uses `0:200` in tests).

### `vol_tree_build` parameters (detailed)

Purpose: build the actual adaptive tree structure and sample source values on box grids.

- `ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta, fun, nd, dpars, zpars, ipars`:
  same meaning as in `vol_tree_mem`; these must match the sizing call.
- `rintl(0:nlevels)` (`real*8`, in): level-wise norm estimates from `vol_tree_mem`.
- `nboxes` (`integer`, in): box count from `vol_tree_mem`.
- `nlevels` (`integer`, in): level count from `vol_tree_mem`.
- `ltree` (`integer`, in): tree array length from `vol_tree_mem`.
- `itree(ltree)` (`integer`, out): packed tree data block.
- `iptr(8)` (`integer`, out): offsets into `itree`.
  - `iptr(1)`: `laddr`
  - `iptr(2)`: `ilevel`
  - `iptr(3)`: `iparent`
  - `iptr(4)`: `nchild`
  - `iptr(5)`: `ichild`
  - `iptr(6)`: `ncoll`
  - `iptr(7)`: `coll`
  - `iptr(8)`: end/total tree length marker
- `centers(ndim,nboxes)` (`real*8`, out): box centers.
- `boxsize(0:nlevels)` (`real*8`, out): physical box size by level.
- `fvals(nd,norder**ndim,nboxes)` (`real*8`, out): sampled source values on tensor-product nodes per box.

### Required call contract

- Always call `vol_tree_mem` first.
- Allocate outputs for `vol_tree_build` using returned `nboxes`, `nlevels`, and `ltree`.
- Keep control and callback parameters consistent between the two calls.
