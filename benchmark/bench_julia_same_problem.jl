using BoxDMK

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

function run_once()
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

    rsign = (rsig * 1.0)^(ndim[] / 2)
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

    t0 = time_ns()
    vol_tree_mem!(ndim, ipoly, iperiod, epstree, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, ifnewtree, nboxes, nlevels, ltree, rintl)

    npbox = Int(norder[])^Int(ndim[])
    itree = Vector{Cint}(undef, Int(ltree[]))
    iptr = Vector{Cint}(undef, 8)
    centers = Vector{Float64}(undef, Int(ndim[] * nboxes[]))
    boxsize = Vector{Float64}(undef, Int(nlevels[] + 1))
    fvals = Vector{Float64}(undef, Int(nd[] * npbox * nboxes[]))

    vol_tree_build!(ndim, ipoly, iperiod, epstree, zk, boxlen, norder, iptype, eta,
        cb, nd, dpars, zpars, ipars, rintl, nboxes, nlevels, ltree,
        itree, iptr, centers, boxsize, fvals)
    t1 = time_ns()

    ntarg = Ref{Cint}(1)
    pot = zeros(Float64, max(1, Int(nd[] * npbox * nboxes[])))
    grad = zeros(Float64, 1)
    hess = zeros(Float64, 1)
    targs = zeros(Float64, Int(ndim[] * ntarg[]))
    pote = zeros(Float64, Int(nd[] * ntarg[]))
    grade = zeros(Float64, max(1, Int(nd[] * ndim[] * ntarg[])))
    hesse = zeros(Float64, 1)
    tottimeinfo = zeros(Float64, 20)

    targs .= 0.0

    bdmk!(nd, ndim, eps, ikernel, beta, ipoly, norder, Ref{Cint}(npbox),
        nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, ifpgh,
        pot, grad, hess, ntarg, targs, ifpghtarg, pote, grade, hesse, tottimeinfo)
    t2 = time_ns()

    pnorm = sqrt(sum(abs2, pot))

    tree_build_s = (t1 - t0) / 1e9
    solve_s = (t2 - t1) / 1e9
    total_s = (t2 - t0) / 1e9

    return (; tree_build_s, solve_s, total_s, eps=eps[], pnorm, nboxes=Int(nboxes[]), nlevels=Int(nlevels[]))
end

function main()
    r = run_once()
    println("BENCH_JULIA tree_build_s=$(round(r.tree_build_s, digits=6)) solve_s=$(round(r.solve_s, digits=6)) total_s=$(round(r.total_s, digits=6)) eps=$(r.eps) pnorm=$(r.pnorm) nboxes=$(r.nboxes) nlevels=$(r.nlevels)")
end

main()
