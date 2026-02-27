using BoxDMK

function density(x, problem)
    r2 = sum(abs2, x)
    return exp(-10r2)
end

function main()

    prob = BDMKProblem(density=density, nd=1, ndim=3, ikernel=1)
    opts = BDMKOptions(eps=1e-6, norder=8)

    @info "Running BoxDMK..."
    tree, res = solve_problem(prob; compute=:potential, opts=opts)

    @info "Evaluating targets..."
    targets = rand(3, 20) .- 0.5
    tres = evaluate_targets(prob, tree, targets; compute=:potential, eps=opts.eps)

    nothing
end

main()
