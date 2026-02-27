module BoxDMK

export set_libboxdmk_path!, libboxdmk_path
export vol_tree_mem!, vol_tree_build!, bdmk!

const _libboxdmk_override = Ref{Union{Nothing,String}}(nothing)

function set_libboxdmk_path!(path::AbstractString)
    _libboxdmk_override[] = String(path)
    return _libboxdmk_override[]
end

function _default_libboxdmk_path()
    return normpath(joinpath(@__DIR__, "..", "..", "build", "libboxdmk.so"))
end

function libboxdmk_path()
    return something(_libboxdmk_override[], _default_libboxdmk_path())
end

function vol_tree_mem!(
    ndim::Ref{Cint}, ipoly::Ref{Cint}, iperiod::Ref{Cint},
    eps::Ref{Cdouble}, zk::Ref{ComplexF64}, boxlen::Ref{Cdouble},
    norder::Ref{Cint}, iptype::Ref{Cint}, eta::Ref{Cdouble},
    funptr::Ptr{Cvoid}, nd::Ref{Cint}, dpars::Vector{Float64},
    zpars::Vector{ComplexF64}, ipars::Vector{Cint}, ifnewtree::Ref{Cint},
    nboxes::Ref{Cint}, nlevels::Ref{Cint}, ltree::Ref{Cint},
    rintl::Vector{Float64},
)
    length(rintl) >= 201 || throw(ArgumentError("rintl must have length >= 201"))

    ccall((:boxdmk_vol_tree_mem, libboxdmk_path()), Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ref{ComplexF64}, Ref{Cdouble},
         Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{Cvoid}, Ref{Cint}, Ptr{Cdouble},
         Ptr{ComplexF64}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{Cdouble}),
        ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta,
        funptr, nd, dpars, zpars, ipars, ifnewtree, nboxes, nlevels, ltree, rintl,
    )
    return nothing
end

function vol_tree_build!(
    ndim::Ref{Cint}, ipoly::Ref{Cint}, iperiod::Ref{Cint},
    eps::Ref{Cdouble}, zk::Ref{ComplexF64}, boxlen::Ref{Cdouble},
    norder::Ref{Cint}, iptype::Ref{Cint}, eta::Ref{Cdouble},
    funptr::Ptr{Cvoid}, nd::Ref{Cint}, dpars::Vector{Float64},
    zpars::Vector{ComplexF64}, ipars::Vector{Cint}, rintl::Vector{Float64},
    nboxes::Ref{Cint}, nlevels::Ref{Cint}, ltree::Ref{Cint},
    itree::Vector{Cint}, iptr::Vector{Cint}, centers::Vector{Float64},
    boxsize::Vector{Float64}, fvals::Vector{Float64},
)
    length(iptr) >= 8 || throw(ArgumentError("iptr must have length >= 8"))

    ccall((:boxdmk_vol_tree_build, libboxdmk_path()), Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ref{ComplexF64}, Ref{Cdouble},
         Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ptr{Cvoid}, Ref{Cint}, Ptr{Cdouble},
         Ptr{ComplexF64}, Ptr{Cint}, Ptr{Cdouble}, Ref{Cint}, Ref{Cint}, Ref{Cint},
         Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        ndim, ipoly, iperiod, eps, zk, boxlen, norder, iptype, eta,
        funptr, nd, dpars, zpars, ipars, rintl, nboxes, nlevels, ltree,
        itree, iptr, centers, boxsize, fvals,
    )
    return nothing
end

function bdmk!(
    nd::Ref{Cint}, ndim::Ref{Cint}, eps::Ref{Cdouble}, ikernel::Ref{Cint},
    beta::Ref{Cdouble}, ipoly::Ref{Cint}, norder::Ref{Cint}, npbox::Ref{Cint},
    nboxes::Ref{Cint}, nlevels::Ref{Cint}, ltree::Ref{Cint},
    itree::Vector{Cint}, iptr::Vector{Cint}, centers::Vector{Float64},
    boxsize::Vector{Float64}, fvals::Vector{Float64}, ifpgh::Ref{Cint},
    pot::Vector{Float64}, grad::Vector{Float64}, hess::Vector{Float64},
    ntarg::Ref{Cint}, targs::Vector{Float64}, ifpghtarg::Ref{Cint},
    pote::Vector{Float64}, grade::Vector{Float64}, hesse::Vector{Float64},
    tottimeinfo::Vector{Float64},
)
    length(iptr) >= 8 || throw(ArgumentError("iptr must have length >= 8"))

    ccall((:boxdmk_bdmk, libboxdmk_path()), Cvoid,
        (Ref{Cint}, Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ref{Cdouble}, Ref{Cint},
         Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cint},
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
         Ptr{Cdouble}, Ref{Cint}, Ptr{Cdouble}, Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble},
         Ptr{Cdouble}, Ptr{Cdouble}),
        nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ltree,
        itree, iptr, centers, boxsize, fvals, ifpgh, pot, grad, hess,
        ntarg, targs, ifpghtarg, pote, grade, hesse, tottimeinfo,
    )
    return nothing
end

end
