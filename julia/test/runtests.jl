using BoxDMK
using Test
using Random
using LinearAlgebra

function rhs_cb_simple(nd::Cint, xyz::Ptr{Cdouble}, dpars::Ptr{Cdouble}, zpars::Ptr{ComplexF64}, ipars::Ptr{Cint}, f::Ptr{Cdouble})::Cvoid
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

# Port of the active Fortran rhs path in test/bdmk/testbdmk.f for ikernel=1, ndim=3.
function rhs_cb_fortran_workflow(nd::Cint, xyz::Ptr{Cdouble}, dpars::Ptr{Cdouble}, zpars::Ptr{ComplexF64}, ipars::Ptr{Cint}, f::Ptr{Cdouble})::Cvoid
    ndim = Int(unsafe_load(ipars, 1))
    ng = Int(unsafe_load(ipars, 3))

    @inbounds for ind in 1:Int(nd)
        fi = 0.0
        for i in 1:ng
            idp = (i - 1) * 5
            rr = 0.0
            for k in 1:ndim
                xk = unsafe_load(xyz, k) - unsafe_load(dpars, idp + k)
                rr += xk * xk
            end
            sigma = unsafe_load(dpars, idp + 4)
            strength = unsafe_load(dpars, idp + 5)
            fi += strength * exp(-rr / sigma) * (-2 * ndim + 4 * rr / sigma) / sigma
        end
        unsafe_store!(f, fi, ind)
    end
    return
end

# Julia equivalent of the uexact path used in test/bdmk/testbdmk.f for ndim=3, ikernel=1.
function uexact_laplace_3d(targ::NTuple{3,Float64}, dpars::Vector{Float64}, ng::Int)
    s = 0.0
    @inbounds for i in 1:ng
        idp = (i - 1) * 5
        dx = targ[1] - dpars[idp + 1]
        dy = targ[2] - dpars[idp + 2]
        dz = targ[3] - dpars[idp + 3]
        sigma = dpars[idp + 4]
        strength = dpars[idp + 5]
        r2 = dx * dx + dy * dy + dz * dz
        s += strength * exp(-r2 / sigma)
    end
    return -4 * pi * s
end

function gauss_legendre_nodes(n::Int)
    if n == 1
        return [0.0]
    end
    b = [k / sqrt(4 * k * k - 1) for k in 1:(n - 1)]
    vals = eigvals(SymTridiagonal(zeros(n), b))
    return sort(collect(real(vals)))
end

function tensor_nodes_3d(n::Int)
    x1d = gauss_legendre_nodes(n)
    np = n^3
    x = zeros(Float64, 3, np)
    ipt = 0
    for i in 1:n
        for j in 1:n
            for k in 1:n
                ipt += 1
                x[1, ipt] = x1d[k]
                x[2, ipt] = x1d[j]
                x[3, ipt] = x1d[i]
            end
        end
    end
    return x
end

@testset "smoke" begin
    @test isdefined(BoxDMK, :vol_tree_mem!)
    @test isdefined(BoxDMK, :vol_tree_build!)
    @test isdefined(BoxDMK, :bdmk!)
    @test isdefined(BoxDMK, :BDMKProblem)
    @test isdefined(BoxDMK, :BDMKOptions)
    @test isdefined(BoxDMK, :BDMKTree)
    @test isdefined(BoxDMK, :BDMKResult)
    @test isdefined(BoxDMK, :YukawaProblem)
    @test isdefined(BoxDMK, :LaplaceProblem)
    @test isdefined(BoxDMK, :SqrtLaplaceProblem)
    @test isdefined(BoxDMK, :build_tree)
    @test isdefined(BoxDMK, :solve)
    @test isdefined(BoxDMK, :evaluate_targets)
    @test isdefined(BoxDMK, :run)
end

@testset "problem-constructors" begin
    density(x, problem) = 1.0

    y = YukawaProblem(density=density, ndim=2, nd=1, beta=2.5)
    l = LaplaceProblem(density=density, ndim=3, nd=1)
    s = SqrtLaplaceProblem(density=density, ndim=2, nd=2, beta=0.75)

    @test y.ikernel == 0
    @test l.ikernel == 1
    @test s.ikernel == 2

    @test y.ipars[2] == Cint(0)
    @test l.ipars[2] == Cint(1)
    @test s.ipars[2] == Cint(2)

    @test y.beta == 2.5
    @test s.beta == 0.75
end

@testset "integration" begin
    libpath = BoxDMK.libboxdmk_path()
    @test isfile(libpath)

    cb = @cfunction(rhs_cb_simple, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}))

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

@testset "fortran-workflow-via-julia-functions" begin
    cb = @cfunction(rhs_cb_fortran_workflow, Cvoid, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}))

    ikernel = Ref{Cint}(1)
    beta = Ref{Cdouble}(1.0)
    eps = Ref{Cdouble}(1e-6)
    rsig = 1e-4

    nd = Ref{Cint}(1)
    ndim = Ref{Cint}(3)
    norder = Ref{Cint}(16)
    ipoly = Ref{Cint}(0)
    ifpgh = Ref{Cint}(1)
    ifpghtarg = Ref{Cint}(0)

    iperiod = Ref{Cint}(0)
    iptype = Ref{Cint}(2)
    eta = Ref{Cdouble}(0.0)
    ifnewtree = Ref{Cint}(0)
    boxlen = Ref{Cdouble}(1.18)
    zk = Ref{ComplexF64}(30.0 + 0.0im)
    epstree = Ref{Cdouble}(eps[] * 500.0)

    ipars = fill(Cint(0), 256)
    dpars = zeros(Float64, 1000)
    zpars = zeros(ComplexF64, 16)

    ipars[1] = ndim[]
    ipars[2] = ikernel[]
    ipars[3] = 2
    ipars[5] = iperiod[]
    ipars[10] = max(ifpgh[], ifpghtarg[])

    delta = 1.0
    rsign = (rsig * delta)^(ndim[] / 2)
    dpars[201] = beta[]

    dpars[1] = 0.1
    dpars[2] = 0.02
    dpars[3] = 0.04
    dpars[4] = rsig
    dpars[5] = 1 / pi / rsign

    dpars[6] = 0.03
    dpars[7] = -0.1
    dpars[8] = 0.05
    dpars[9] = rsig / 2
    dpars[10] = -0.5 / pi / rsign

    nboxes = Ref{Cint}(0)
    nlevels = Ref{Cint}(0)
    ltree = Ref{Cint}(0)
    rintl = zeros(Float64, 201)

    BoxDMK.vol_tree_mem!(ndim, ipoly, iperiod, epstree, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, ifnewtree, nboxes, nlevels, ltree, rintl)

    @test nboxes[] > 100
    @test nlevels[] >= 3
    @test ltree[] > 0

    npbox = Int(norder[])^Int(ndim[])
    itree = Vector{Cint}(undef, Int(ltree[]))
    iptr = Vector{Cint}(undef, 8)
    centers = Vector{Float64}(undef, Int(ndim[] * nboxes[]))
    boxsize = Vector{Float64}(undef, Int(nlevels[] + 1))
    fvals = Vector{Float64}(undef, Int(nd[] * npbox * nboxes[]))

    BoxDMK.vol_tree_build!(ndim, ipoly, iperiod, epstree, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, rintl, nboxes, nlevels, ltree,
        itree, iptr, centers, boxsize, fvals)

    nhess = Int(ndim[] * (ndim[] + 1) ÷ 2)
    ntarg = Ref{Cint}(20)
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

    laddr_ptr = 1
    nchild_ptr = Int(iptr[4])
    nlfbox = 0
    for ilevel in 1:Int(nlevels[])
        first_box = Int(itree[2 * ilevel + 1])
        last_box = Int(itree[2 * ilevel + 2])
        for ibox in first_box:last_box
            if itree[nchild_ptr + ibox - 1] == 0
                nlfbox += 1
            end
        end
    end

    ntot = npbox * nlfbox * Int(ifpgh[]) + Int(ntarg[]) * Int(ifpghtarg[])
    pnorm = sqrt(sum(abs2, pot))
    pref = zeros(Float64, Int(ntarg[]))
    for i in 1:Int(ntarg[])
        targ = (
            targs[3 * (i - 1) + 1],
            targs[3 * (i - 1) + 2],
            targs[3 * (i - 1) + 3],
        )
        pref[i] = uexact_laplace_3d(targ, dpars, Int(ipars[3]))
    end
    centers_m = reshape(centers, Int(ndim[]), Int(nboxes[]))
    pot_m = reshape(pot, Int(nd[]), npbox, Int(nboxes[]))
    xref = tensor_nodes_3d(Int(norder[]))

    abserr2 = 0.0
    rnorm2 = 0.0
    for ilevel in 0:Int(nlevels[])
        first_box = Int(itree[2 * ilevel + 1])
        last_box = Int(itree[2 * ilevel + 2])
        bs = boxsize[ilevel + 1] / 2
        for ibox in first_box:last_box
            if itree[nchild_ptr + ibox - 1] == 0
                for j in 1:npbox
                    targ = (
                        centers_m[1, ibox] + xref[1, j] * bs,
                        centers_m[2, ibox] + xref[2, j] * bs,
                        centers_m[3, ibox] + xref[3, j] * bs,
                    )
                    pex = uexact_laplace_3d(targ, dpars, Int(ipars[3]))
                    pnum = pot_m[1, j, ibox]
                    rnorm2 += pex * pex
                    abserr2 += (pex - pnum) * (pex - pnum)
                end
            end
        end
    end

    abserrp = sqrt(abserr2)
    rnormp = sqrt(rnorm2)
    if rnormp < 1e-6
        perr = abserrp
    else
        perr = abserrp / rnormp
    end
    @test nlfbox > 0
    @test ntot > 0
    @test pnorm > 1e6
    @test all(isfinite, tottimeinfo)
    @test isfinite(perr)
    @test perr < 1e-5
end

@testset "highlevel-api" begin
    rsig = 5e-2
    ndim = 2
    nd = 1
    rsign = rsig^(ndim / 2)

    dpars = zeros(Float64, 1000)
    dpars[1] = 0.1
    dpars[2] = -0.1
    dpars[3] = 0.0
    dpars[4] = rsig
    dpars[5] = 1 / pi / rsign
    dpars[201] = 1.0

    ipars = fill(Cint(0), 256)
    ipars[3] = 1
    ipars[5] = 0

    function density_highlevel(x::Vector{Float64}, problem::BDMKProblem)
        ng = Int(problem.ipars[3])
        fi = 0.0
        for i in 1:ng
            idp = (i - 1) * 5
            rr = 0.0
            for k in 1:problem.ndim
                rr += (x[k] - problem.dpars[idp + k])^2
            end
            sigma = problem.dpars[idp + 4]
            strength = problem.dpars[idp + 5]
            fi += strength * exp(-rr / sigma) * (-2 * problem.ndim + 4 * rr / sigma) / sigma
        end
        return [fi]
    end

    problem = BDMKProblem(
        density=density_highlevel,
        nd=1,
        ndim=ndim,
        ikernel=1,
        beta=1.0,
        boxlen=1.18,
        dpars=dpars,
        zpars=zeros(ComplexF64, 16),
        ipars=ipars,
    )
    opts = BDMKOptions(eps=1e-2, norder=2, ipoly=0, iptype=2, eta=0.0, epstree_factor=500.0)

    tree = build_tree(problem, opts)
    @test tree.nboxes > 0
    @test tree.nlevels >= 0
    @test tree.npbox == opts.norder^problem.ndim

    result = solve(problem, tree; compute=:potential, eps=opts.eps)
    @test size(result.pot) == (1, tree.npbox, tree.nboxes)
    @test result.pote === nothing

    Random.seed!(1234)
    targets = rand(problem.ndim, 4) .- 0.5
    tres = evaluate_targets(problem, tree, targets; compute=:potential, eps=opts.eps)
    @test tres.pote !== nothing
    @test size(tres.pote) == (1, size(targets, 2))

    tree2, runres = BoxDMK.run(problem; targets=targets, compute=:potential, opts=opts)
    @test tree2.nboxes == tree.nboxes
    @test size(runres.pote) == (1, size(targets, 2))

    pot_ll = zeros(Float64, problem.nd * tree.npbox * tree.nboxes)
    grad_ll = zeros(Float64, 1)
    hess_ll = zeros(Float64, 1)
    ntarg = Ref{Cint}(1)
    targs = zeros(Float64, problem.ndim)
    ifpghtarg = Ref{Cint}(0)
    pote_ll = zeros(Float64, 1)
    grade_ll = zeros(Float64, 1)
    hesse_ll = zeros(Float64, 1)
    timeinfo = zeros(Float64, 20)

    bdmk!(
        Ref{Cint}(problem.nd), Ref{Cint}(problem.ndim), Ref{Cdouble}(opts.eps), Ref{Cint}(problem.ikernel),
        Ref{Cdouble}(problem.beta), Ref{Cint}(opts.ipoly), Ref{Cint}(tree.norder), Ref{Cint}(tree.npbox),
        Ref{Cint}(tree.nboxes), Ref{Cint}(tree.nlevels), Ref{Cint}(tree.ltree), tree.itree, tree.iptr,
        tree.centers, tree.boxsize, tree.fvals, Ref{Cint}(1), pot_ll, grad_ll, hess_ll,
        ntarg, targs, ifpghtarg, pote_ll, grade_ll, hesse_ll, timeinfo
    )

    pot_ll_shaped = reshape(pot_ll, 1, tree.npbox, tree.nboxes)
    rel = norm(vec(result.pot .- pot_ll_shaped)) / norm(vec(pot_ll_shaped))
    @test rel < 1e-12
end
