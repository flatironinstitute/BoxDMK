# High-level Julia interface on top of low-level ccall wrappers.

export BDMKProblem, BDMKOptions, BDMKTree, BDMKResult
export YukawaProblem, LaplaceProblem, SqrtLaplaceProblem
export build_tree, solve, evaluate_targets, run

"""
    BDMKProblem(; density, nd=1, ndim=3, ikernel=1, beta=1.0, boxlen=1.18,
                dpars=zeros(Float64,1000), zpars=zeros(ComplexF64,16), ipars=fill(Cint(0),256))

Problem definition for BDMK solves.

Parameters:
- `density`: Julia callback `density(x, problem)` used to fill source values at tree nodes.
- `nd`: number of output components in the field/source value.
- `ndim`: spatial dimension (`1`, `2`, or `3`).
- `ikernel`: kernel selector (`0`: Yukawa, `1`: Laplace, `2`: square-root Laplace).
- `beta`: kernel parameter used by low-level solver.
- `boxlen`: side length of the computational box.
- `dpars`: real parameter buffer passed to low-level code/callback.
- `zpars`: complex parameter buffer passed to low-level code/callback.
- `ipars`: integer parameter buffer passed to low-level code/callback.
"""
struct BDMKProblem
    nd::Int
    ndim::Int
    ikernel::Int
    beta::Float64
    boxlen::Float64
    density::Function
    dpars::Vector{Float64}
    zpars::Vector{ComplexF64}
    ipars::Vector{Cint}
end

function BDMKProblem(; density::Function,
    nd::Integer=1,
    ndim::Integer=3,
    ikernel::Integer=1,
    beta::Real=1.0,
    boxlen::Real=1.18,
    dpars::Vector{Float64}=zeros(Float64, 1000),
    zpars::Vector{ComplexF64}=zeros(ComplexF64, 16),
    ipars::Vector{Cint}=fill(Cint(0), 256))

    nd > 0 || throw(ArgumentError("nd must be positive"))
    ndim in (1, 2, 3) || throw(ArgumentError("ndim must be 1, 2, or 3"))

    ipars2 = copy(ipars)
    length(ipars2) >= 10 || resize!(ipars2, 10)
    ipars2[1] = Cint(ndim)
    ipars2[2] = Cint(ikernel)

    return BDMKProblem(Int(nd), Int(ndim), Int(ikernel), Float64(beta), Float64(boxlen),
        density, copy(dpars), copy(zpars), ipars2)
end

"""
    YukawaProblem(; kwargs...)

Shorthand for `BDMKProblem(...; ikernel=0, ...)`.
"""
function YukawaProblem(; kwargs...)
    return BDMKProblem(; ikernel=0, kwargs...)
end

"""
    LaplaceProblem(; kwargs...)

Shorthand for `BDMKProblem(...; ikernel=1, ...)`.
"""
function LaplaceProblem(; kwargs...)
    return BDMKProblem(; ikernel=1, kwargs...)
end

"""
    SqrtLaplaceProblem(; kwargs...)

Shorthand for `BDMKProblem(...; ikernel=2, ...)`.
"""
function SqrtLaplaceProblem(; kwargs...)
    return BDMKProblem(; ikernel=2, kwargs...)
end

"""
    BDMKOptions(; eps=1e-6, norder=16, ipoly=0, iptype=2, eta=0.0,
                epstree_factor=500.0, ifnewtree=0, iperiod=0, zk=30+0im)

Options controlling tree construction and solver evaluation.

Parameters:
- `eps`: solve tolerance used in `bdmk!`.
- `norder`: interpolation order per dimension (`npbox = norder^ndim`).
- `ipoly`: interpolation polynomial flag for low-level tree routines.
- `iptype`: point-type flag for low-level tree routines.
- `eta`: tuning parameter forwarded to low-level tree routines.
- `epstree_factor`: scale factor for tree build tolerance (`eps_tree = eps * epstree_factor`).
- `ifnewtree`: low-level tree build mode flag.
- `iperiod`: periodicity flag.
- `zk`: complex kernel parameter forwarded to tree routines.
"""
Base.@kwdef struct BDMKOptions
    eps::Float64 = 1e-6
    norder::Int = 16
    ipoly::Int = 0
    iptype::Int = 2
    eta::Float64 = 0.0
    epstree_factor::Float64 = 500.0
    ifnewtree::Int = 0
    iperiod::Int = 0
    zk::ComplexF64 = 30.0 + 0.0im
end

"""
Container holding the built tree and source-node arrays used by `bdmk!`.
"""
struct BDMKTree
    nd::Int
    ndim::Int
    norder::Int
    npbox::Int
    nboxes::Int
    nlevels::Int
    ltree::Int
    itree::Vector{Cint}
    iptr::Vector{Cint}
    centers::Vector{Float64}
    boxsize::Vector{Float64}
    fvals::Vector{Float64}
    rintl::Vector{Float64}
end

"""
Container of node/target outputs returned by `solve` and `evaluate_targets`.
"""
struct BDMKResult
    pot
    grad
    hess
    pote
    grade
    hesse
    tottimeinfo::Vector{Float64}
    ifpgh::Int
    ifpghtarg::Int
end

const _HL_CALLBACK_CONTEXT = Ref{Any}(nothing)

function _hl_density_trampoline(nd::Cint, xyz::Ptr{Cdouble}, dpars::Ptr{Cdouble}, zpars::Ptr{ComplexF64}, ipars::Ptr{Cint}, fout::Ptr{Cdouble})::Cvoid
    ctx = _HL_CALLBACK_CONTEXT[]
    if ctx === nothing
        @inbounds for i in 1:Int(nd)
            unsafe_store!(fout, 0.0, i)
        end
        return
    end

    problem = ctx::BDMKProblem
    x = Vector{Float64}(undef, problem.ndim)
    @inbounds for k in 1:problem.ndim
        x[k] = unsafe_load(xyz, k)
    end

    vals = nothing
    try
        vals = problem.density(x, problem)
    catch
        vals = 0.0
    end

    if vals isa Number
        v = Float64(vals)
        @inbounds for i in 1:Int(nd)
            unsafe_store!(fout, v, i)
        end
    else
        nout = min(Int(nd), length(vals))
        @inbounds for i in 1:nout
            unsafe_store!(fout, Float64(vals[i]), i)
        end
        @inbounds for i in (nout + 1):Int(nd)
            unsafe_store!(fout, 0.0, i)
        end
    end
    return
end

const _HL_DENSITY_CFUN = Ref{Ptr{Cvoid}}(C_NULL)

function __init__()
    _HL_DENSITY_CFUN[] = @cfunction(_hl_density_trampoline, Cvoid,
        (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}))
end

function _compute_flag(compute)
    c = Symbol(compute)
    if c == :potential
        return 1
    elseif c == :gradient
        return 2
    elseif c == :hessian
        return 3
    else
        throw(ArgumentError("compute must be :potential, :gradient, or :hessian"))
    end
end

function _prepare_problem_ipars(problem::BDMKProblem, ifout::Int, iperiod::Int)
    ipars = copy(problem.ipars)
    length(ipars) >= 10 || resize!(ipars, 10)
    ipars[1] = Cint(problem.ndim)
    ipars[2] = Cint(problem.ikernel)
    ipars[5] = Cint(iperiod)
    ipars[10] = Cint(ifout)
    return ipars
end

function build_tree(problem::BDMKProblem, opts::BDMKOptions=BDMKOptions())
    nd = Ref{Cint}(problem.nd)
    ndim = Ref{Cint}(problem.ndim)
    ipoly = Ref{Cint}(opts.ipoly)
    iperiod = Ref{Cint}(opts.iperiod)
    norder = Ref{Cint}(opts.norder)
    iptype = Ref{Cint}(opts.iptype)
    ifnewtree = Ref{Cint}(opts.ifnewtree)

    eps_tree = Ref{Cdouble}(opts.eps * opts.epstree_factor)
    eta = Ref{Cdouble}(opts.eta)
    zk = Ref{ComplexF64}(opts.zk)
    boxlen = Ref{Cdouble}(problem.boxlen)

    nboxes = Ref{Cint}(0)
    nlevels = Ref{Cint}(0)
    ltree = Ref{Cint}(0)
    rintl = zeros(Float64, 201)

    ipars = _prepare_problem_ipars(problem, 1, opts.iperiod)

    _HL_CALLBACK_CONTEXT[] = problem
    vol_tree_mem!(ndim, ipoly, iperiod, eps_tree, zk, boxlen, norder, iptype, eta,
        _HL_DENSITY_CFUN[], nd, problem.dpars, problem.zpars, ipars, ifnewtree,
        nboxes, nlevels, ltree, rintl)

    npbox = Int(opts.norder)^Int(problem.ndim)
    itree = Vector{Cint}(undef, Int(ltree[]))
    iptr = Vector{Cint}(undef, 8)
    centers = Vector{Float64}(undef, Int(problem.ndim * nboxes[]))
    boxsize = Vector{Float64}(undef, Int(nlevels[] + 1))
    fvals = Vector{Float64}(undef, Int(problem.nd * npbox * nboxes[]))

    _HL_CALLBACK_CONTEXT[] = problem
    vol_tree_build!(ndim, ipoly, iperiod, eps_tree, zk, boxlen, norder, iptype, eta,
        _HL_DENSITY_CFUN[], nd, problem.dpars, problem.zpars, ipars, rintl,
        nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals)

    return BDMKTree(problem.nd, problem.ndim, opts.norder, npbox, Int(nboxes[]), Int(nlevels[]),
        Int(ltree[]), itree, iptr, centers, boxsize, fvals, rintl)
end

function solve(problem::BDMKProblem, tree::BDMKTree; compute=:potential, targets=nothing, eps::Real=1e-6)
    ifpgh = _compute_flag(compute)

    nd = Ref{Cint}(problem.nd)
    ndim = Ref{Cint}(problem.ndim)
    eps_ref = Ref{Cdouble}(Float64(eps))
    ikernel = Ref{Cint}(problem.ikernel)
    beta = Ref{Cdouble}(problem.beta)
    ipoly = Ref{Cint}(0)
    norder = Ref{Cint}(tree.norder)
    npbox = Ref{Cint}(tree.npbox)
    nboxes = Ref{Cint}(tree.nboxes)
    nlevels = Ref{Cint}(tree.nlevels)
    ltree = Ref{Cint}(tree.ltree)
    ifpgh_ref = Ref{Cint}(ifpgh)

    nhess = Int(problem.ndim * (problem.ndim + 1) ÷ 2)

    pot_vec = zeros(Float64, problem.nd * tree.npbox * tree.nboxes)
    grad_vec = if ifpgh >= 2
        zeros(Float64, problem.nd * problem.ndim * tree.npbox * tree.nboxes)
    else
        zeros(Float64, 1)
    end
    hess_vec = if ifpgh >= 3
        zeros(Float64, problem.nd * nhess * tree.npbox * tree.nboxes)
    else
        zeros(Float64, 1)
    end

    if targets === nothing
        ntarg = Ref{Cint}(1)
        targs_vec = zeros(Float64, max(1, problem.ndim))
        ifpghtarg = Ref{Cint}(0)
        pote_vec = zeros(Float64, max(1, problem.nd))
        grade_vec = zeros(Float64, 1)
        hesse_vec = zeros(Float64, 1)
    else
        size(targets, 1) == problem.ndim || throw(ArgumentError("targets must have size (ndim, ntarg)"))
        nt = size(targets, 2)
        ntarg = Ref{Cint}(nt)
        targs_vec = vec(copy(targets))
        ifpghtarg = Ref{Cint}(ifpgh)
        pote_vec = zeros(Float64, problem.nd * nt)
        grade_vec = if ifpgh >= 2
            zeros(Float64, problem.nd * problem.ndim * nt)
        else
            zeros(Float64, 1)
        end
        hesse_vec = if ifpgh >= 3
            zeros(Float64, problem.nd * nhess * nt)
        else
            zeros(Float64, 1)
        end
    end

    tottimeinfo = zeros(Float64, 20)

    bdmk!(nd, ndim, eps_ref, ikernel, beta, ipoly, norder, npbox,
        nboxes, nlevels, ltree, tree.itree, tree.iptr, tree.centers,
        tree.boxsize, tree.fvals, ifpgh_ref, pot_vec, grad_vec, hess_vec,
        ntarg, targs_vec, ifpghtarg, pote_vec, grade_vec, hesse_vec, tottimeinfo)

    pot = reshape(pot_vec, problem.nd, tree.npbox, tree.nboxes)
    grad = if ifpgh >= 2
        reshape(grad_vec, problem.nd, problem.ndim, tree.npbox, tree.nboxes)
    else
        nothing
    end
    hess = if ifpgh >= 3
        reshape(hess_vec, problem.nd, nhess, tree.npbox, tree.nboxes)
    else
        nothing
    end

    pote = if ifpghtarg[] >= 1
        reshape(pote_vec, problem.nd, Int(ntarg[]))
    else
        nothing
    end
    grade = if ifpghtarg[] >= 2
        reshape(grade_vec, problem.nd, problem.ndim, Int(ntarg[]))
    else
        nothing
    end
    hesse = if ifpghtarg[] >= 3
        reshape(hesse_vec, problem.nd, nhess, Int(ntarg[]))
    else
        nothing
    end

    return BDMKResult(pot, grad, hess, pote, grade, hesse, tottimeinfo, ifpgh, Int(ifpghtarg[]))
end

function evaluate_targets(problem::BDMKProblem, tree::BDMKTree, targets; compute=:potential, eps::Real=1e-6)
    return solve(problem, tree; compute=compute, targets=targets, eps=eps)
end

function run(problem::BDMKProblem; targets=nothing, compute=:potential, opts=BDMKOptions())
    tree = build_tree(problem, opts)
    result = solve(problem, tree; compute=compute, targets=targets, eps=opts.eps)
    return tree, result
end
