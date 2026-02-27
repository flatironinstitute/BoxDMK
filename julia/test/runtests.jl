using Test

include("../src/BoxDMK.jl")
using .BoxDMK

function rhs_cb(nd::Cint, xyz::Ptr{Cdouble}, dpars::Ptr{Cdouble}, zpars::Ptr{ComplexF64}, ipars::Ptr{Cint}, f::Ptr{Cdouble})::Cvoid
    ndim = unsafe_load(ipars, 1)
    r2 = 0.0
    @inbounds for k in 1:ndim
        xk = unsafe_load(xyz, k)
        r2 += xk * xk
    end
    val = exp(-10.0 * r2)
    @inbounds for i in 1:nd
        unsafe_store!(f, val, i)
    end
    return
end

@testset "smoke" begin
    @test isdefined(BoxDMK, :vol_tree_mem!)
    @test isdefined(BoxDMK, :vol_tree_build!)
    @test isdefined(BoxDMK, :bdmk!)
end

@testset "integration" begin
    libpath = BoxDMK.libboxdmk_path()
    @test isfile(libpath)

    cb = @cfunction(rhs_cb, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}))

    ndim = Ref{Cint}(2)
    ipoly = Ref{Cint}(0)
    iperiod = Ref{Cint}(0)
    eps = Ref{Cdouble}(1e-1)
    zk = Ref{ComplexF64}(30.0 + 0.0im)
    boxlen = Ref{Cdouble}(1.18)
    norder = Ref{Cint}(3)
    iptype = Ref{Cint}(2)
    eta = Ref{Cdouble}(0.0)
    nd = Ref{Cint}(1)
    ifnewtree = Ref{Cint}(0)

    ipars = fill(Cint(0), 32)
    ipars[1] = ndim[]
    dpars = zeros(Float64, 32)
    zpars = zeros(ComplexF64, 8)

    nboxes = Ref{Cint}(0)
    nlevels = Ref{Cint}(0)
    ltree = Ref{Cint}(0)
    rintl = zeros(Float64, 201)

    BoxDMK.vol_tree_mem!(ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, ifnewtree, nboxes, nlevels, ltree, rintl)

    @test nboxes[] > 0
    @test nlevels[] >= 0
    @test ltree[] > 0

    npbox = Int(norder[])^Int(ndim[])
    itree = Vector{Cint}(undef, Int(ltree[]))
    iptr = Vector{Cint}(undef, 8)
    centers = Vector{Float64}(undef, Int(ndim[] * nboxes[]))
    boxsize = Vector{Float64}(undef, Int(nlevels[] + 1))
    fvals = Vector{Float64}(undef, Int(nd[] * npbox * nboxes[]))

    BoxDMK.vol_tree_build!(ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, rintl, nboxes, nlevels, ltree,
        itree, iptr, centers, boxsize, fvals)

    @test all(isfinite, fvals)
    @test sum(abs, fvals) > 0.0

    ikernel = Ref{Cint}(1)
    beta = Ref{Cdouble}(1.0)
    ifpgh = Ref{Cint}(1)
    ntarg = Ref{Cint}(1)
    ifpghtarg = Ref{Cint}(0)

    pot = zeros(Float64, Int(nd[] * npbox * nboxes[]))
    grad = zeros(Float64, 1)
    hess = zeros(Float64, 1)
    targs = zeros(Float64, Int(ndim[] * ntarg[]))
    pote = zeros(Float64, max(1, Int(nd[] * ntarg[])))
    grade = zeros(Float64, 1)
    hesse = zeros(Float64, 1)
    tottimeinfo = zeros(Float64, 20)

    BoxDMK.bdmk!(nd, ndim, eps, ikernel, beta, ipoly, norder, Ref{Cint}(npbox),
        nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, ifpgh,
        pot, grad, hess, ntarg, targs, ifpghtarg, pote, grade, hesse, tottimeinfo)

    @test all(isfinite, pot)
    @test sum(abs, pot) > 0.0
end
