# Box DMK

Box DMK is a tool for evaluating volume integrals of Yukawa, Laplace, and square-root Laplace kernels. It is based on the Fortran code developed by Shidong, and provides a Julia wrapper created by Xuanzhao. Original paper of DMK: https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.22240.

## Documentation

- Julia low-level interface: `docs/julia-interface.md`

## Build and Test (Make)

The `all` target compiles and runs the test program (`./int2-bdmk`).

```bash
make clean
make HOST=linux-gfortran FEND='-lopenblas' all
```

Notes:
- Default `HOST=linux-ifort` expects `ifx` and Intel MKL.
- If MKL is installed, you can use:

```bash
make clean
make HOST=linux-gfortran all
```

## Build and Test (CMake)

Load dependencies with environment modules:

```bash
module purge
module load gcc cmake lib/openblas
```

Configure and build (default: OpenBLAS):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Optional MKL support (if available in your environment):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBOXDMK_USE_MKL=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Julia Interface

Build the shared library:

```bash
module purge
module load gcc cmake lib/openblas
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Run Julia smoke tests:

```bash
julia --project=julia julia/test/runtests.jl
```
