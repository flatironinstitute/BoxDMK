# Julia CCall Interface Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a low-level Julia interface that calls `vol_tree_mem`, `vol_tree_build`, and `bdmk` through a stable Fortran `bind(C)` ABI.

**Architecture:** Add a Fortran shim layer exporting C ABI symbols, compile a shared library via CMake, and expose thin Julia wrappers using `ccall`. Keep the API low-level and explicit with caller-managed buffers and callback pointers.

**Tech Stack:** Fortran (iso_c_binding), CMake, Julia (`ccall`, `Test` stdlib)

---

### Task 1: Add failing Julia smoke tests for planned low-level API

**Files:**
- Create: `julia/Project.toml`
- Create: `julia/test/runtests.jl`
- Test: `julia/test/runtests.jl`

**Step 1: Write the failing test**

```julia
using Test
include("../src/BoxDMK.jl")
using .BoxDMK

@test isdefined(BoxDMK, :vol_tree_mem!)
@test isdefined(BoxDMK, :vol_tree_build!)
@test isdefined(BoxDMK, :bdmk!)
```

**Step 2: Run test to verify it fails**

Run: `julia --project=julia julia/test/runtests.jl`
Expected: FAIL because `julia/src/BoxDMK.jl` does not exist yet.

**Step 3: Commit**

```bash
git add julia/Project.toml julia/test/runtests.jl
git commit -m "test: add failing Julia smoke tests for low-level wrappers"
```

### Task 2: Add Fortran `bind(C)` shim exports for the 3 routines

**Files:**
- Create: `src/bdmk/bdmk_c_api.f90`
- Modify: `CMakeLists.txt`
- Test: shared library build

**Step 1: Write the failing test command**

Run: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`
Expected: FAIL for unresolved `boxdmk_*` symbols (wrappers not yet implemented/linked).

**Step 2: Write minimal implementation**

- Implement `boxdmk_vol_tree_mem`, `boxdmk_vol_tree_build`, `boxdmk_bdmk` in `src/bdmk/bdmk_c_api.f90`.
- Add callback trampoline for `fun` pointer using `iso_c_binding`.
- Update CMake to build `boxdmk` shared library containing shim + core sources.

**Step 3: Run test to verify it passes**

Run: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`
Expected: PASS and `build/libboxdmk.so` exists.

**Step 4: Commit**

```bash
git add src/bdmk/bdmk_c_api.f90 CMakeLists.txt
git commit -m "feat: export C ABI shims for core bdmk routines"
```

### Task 3: Implement Julia low-level wrappers using `ccall`

**Files:**
- Create: `julia/src/BoxDMK.jl`
- Modify: `julia/test/runtests.jl`
- Test: `julia/test/runtests.jl`

**Step 1: Run failing tests**

Run: `julia --project=julia julia/test/runtests.jl`
Expected: FAIL because wrappers or module implementation are incomplete.

**Step 2: Write minimal implementation**

- Implement `module BoxDMK`.
- Add `libboxdmk_path()` and `set_libboxdmk_path!`.
- Add `vol_tree_mem!`, `vol_tree_build!`, `bdmk!` wrappers mapped 1:1 to exported C ABI symbols.
- Keep arguments explicit (`Ref`, `Vector`, `Ptr`) with no hidden allocations.

**Step 3: Run tests to verify pass**

Run: `julia --project=julia julia/test/runtests.jl`
Expected: PASS for symbol/wrapper smoke tests.

**Step 4: Commit**

```bash
git add julia/src/BoxDMK.jl julia/test/runtests.jl
git commit -m "feat: add Julia low-level ccall wrappers for boxdmk"
```

### Task 4: Document usage and verify end-to-end checks

**Files:**
- Modify: `docs/user-interface.md`
- Modify: `README.md`
- Test: CMake build + Julia tests

**Step 1: Write docs updates**

- Add Julia section describing low-level API scope.
- Document callback contract and required buffer shapes.
- Link new Julia interface doc path from `README.md`.

**Step 2: Run verification**

Run:
`module purge && module load gcc cmake lib/openblas && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build --output-on-failure && julia --project=julia julia/test/runtests.jl`
Expected: build passes, ctest passes, Julia smoke tests pass.

**Step 3: Commit**

```bash
git add README.md docs/user-interface.md
git commit -m "docs: add Julia low-level interface usage notes"
```
