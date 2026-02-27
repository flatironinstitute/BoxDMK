# Julia Low-Level Interface Design (boxdmk)

## Goal
Provide a first-pass Julia interface that imports and calls the low-level Fortran routines `vol_tree_mem`, `vol_tree_build`, and `bdmk` via `ccall`.

## Scope
- Add a stable C ABI surface in Fortran using `bind(C)` shims.
- Build a shared library (`libboxdmk.so`) that exports only the shim symbols needed by Julia.
- Add a Julia module exposing low-level wrappers with explicit buffer arguments.
- Document callback and array-shape requirements.

Out of scope for this first pass:
- High-level Julia convenience API.
- Automatic allocation of tree/solver buffers.
- End-to-end scientific validation from Julia.

## Architecture

### 1) Fortran C-ABI shim layer
A new shim source file provides `bind(C)` entry points:
- `boxdmk_vol_tree_mem`
- `boxdmk_vol_tree_build`
- `boxdmk_bdmk`

The shim translates C-interoperable scalar/array types to existing internal routine signatures.

For `vol_tree_mem` / `vol_tree_build`, the user callback is accepted as a C function pointer and invoked through a Fortran trampoline.

### 2) Build system updates
CMake adds a shared library target containing:
- existing core Fortran sources
- new C-ABI shim source

This produces `libboxdmk.so` for Julia `ccall`.

### 3) Julia low-level module
A new Julia module:
- loads `libboxdmk` from a configurable path
- provides thin wrappers around the 3 exported symbols
- uses explicit `Ref`/`Vector` arguments for all inputs and outputs

The API remains intentionally low-level and mirrors Fortran calling style.

## Data Flow
1. Julia user prepares callback and parameter buffers.
2. Julia calls `vol_tree_mem` wrapper to get `nboxes`, `nlevels`, `ltree`, `rintl`.
3. Julia allocates tree arrays (`itree`, `iptr`, `centers`, `boxsize`, `fvals`).
4. Julia calls `vol_tree_build` wrapper to build the adaptive tree.
5. Julia allocates output arrays and calls `bdmk` wrapper.

## Error Handling
- This first pass keeps legacy behavior and does not change Fortran error model.
- Wrapper enforces no extra policy; it assumes caller provides correctly sized buffers.
- Documentation calls out minimum required buffer sizes.

## Testing Strategy
- Build-level verification: shared library target builds successfully.
- Julia smoke tests:
  - module loads
  - expected symbols/wrapper functions exist
  - callback type and call signatures compile

## Risks and Mitigations
- ABI mismatch risk: mitigated by centralizing exported ABI in `bind(C)` shims.
- Callback interoperability risk: mitigated via explicit `iso_c_binding` trampoline.
- Buffer sizing misuse: mitigated with docs and strict low-level API naming.
