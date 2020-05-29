"""
    ph_subproblem_solve(pb, id_scen, u_scen, x_scen, μ, params)

Solve subproblem associated with scenario `id_scen`, primal point `x_scen`, dual point `u_scen`
and regularization `μ`.
"""
function ph_subproblem_solve(pb::Problem, id_scen::ScenarioId, u_scen, x_scen, μ, params::AbstractDict)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(optimizer_with_attributes(params[:optimizer], params[:optimizer_params]...))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)

    # Augmented lagragian subproblem full objective
    obj += dot(u_scen, y) + (1/(2*μ)) * sum((y[i]-x_scen[i])^2 for i in 1:n)
    @objective(model, MOI.MIN_SENSE, obj)

    optimize!(model)
    if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
        @warn "Subproblem status: (scenario $id_scen) " primal_status(model) dual_status(model) termination_status(model)
    end

    return JuMP.value.(y)
end


"""
    ph_initialization!(x, u, y, pb, μ, subpbparams, printlev)

Compute a first global feasible point by solving each scenario independently and projecting
the global strategy obtained onto the non-anticipatory subspace.
"""
function ph_initialization!(x, u, y, pb, μ, subpbparams, printlev)
    printlev>0 && print("Initialisation... ")

    # Subproblem solves
    for id_scen in 1:pb.nscenarios
        y[id_scen, :] = ph_subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, subpbparams)
    end

    # projection on non anticipatory subspace
    nonanticipatory_projection!(x, pb, y)

    # multiplier update
    u[:] = ((1/μ) * (y-x))[:]

    printlev>0 && println("done")
    return
end

function ph_print_log(pb, x, u, printlev, hist, it, nscenariostreated, primres, dualres, computingtime, tinit, callback)
    objval = objective_value(pb, x)
    dot_xu = dot(pb, x, u)

    printlev>0 && @printf "%3i   %.2e       %.10e  %.10e   % .3e % .16e\n" it nscenariostreated primres dualres dot_xu objval

    if hist!==nothing
        push!(hist[:functionalvalue], objval)
        push!(hist[:computingtime], computingtime)
        push!(hist[:time], time() - tinit)
        push!(hist[:primres], primres)
        push!(hist[:dualres], dualres)
        push!(hist[:residual], sqrt(primres^2+dualres^2))
        haskey(hist, :approxsol) && size(hist[:approxsol])==size(x) && push!(hist[:dist_opt], norm(hist[:approxsol] - x))
    end

    if (callback!==nothing)
        callback(pb, x, hist)
    end
    return
end

function ph_hasconverged(pb, x, primres, dualres, ϵ_abs, ϵ_rel)
    ϵ = ϵ_abs + ϵ_rel * norm(pb, x)
    return (primres < ϵ) && (dualres < ϵ)
end

"""
    x = solve_progressivehedging(pb::Problem)

Run the classical Progressive Hedging scheme on problem `pb`.

Stopping criterion is based on primal dual residual, maximum iterations or time
can also be set. Return a feasible point `x`.

## Keyword arguments:
- `ϵ_abs`: Absolute tolerance on residual.
- `ϵ_rel`: Relative tolerance on residual.
- `μ`: Regularization parameter.
- `maxtime`: Limit on time spent in `solve_progressivehedging`.
- `maxcomputingtime`: Limit time spent in computations, excluding computation of initial feasible point and computations required by plottting / logging.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
- `hist`: if not `nothing`, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{String, Any}` storing parameters for the optimizer.
- `callback`: either nothing or a function `callback(pb, x, hist)::nothing` called at each log phase. `x` is the current feasible global iterate.
"""
function solve_progressivehedging(pb::Problem; ϵ_abs = 1e-8,
                                               ϵ_rel = 1e-4,
                                               μ = 3,
                                               maxtime = 3600,
                                               maxcomputingtime = Inf,
                                               maxiter = 1e3,
                                               printlev = 1,
                                               printstep = 1,
                                               hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                               optimizer = Ipopt.Optimizer,
                                               optimizer_params = Dict{String, Any}("print_level"=>0),
                                               callback = nothing,
                                               kwargs...)
    display_algopb_stats(pb, "Progressive Hedging", printlev, ϵ_abs=ϵ_abs, ϵ_rel=ϵ_rel, μ=μ, maxtime=maxtime, maxcomputingtime=maxcomputingtime, maxiter=maxiter)

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)

    y = zeros(Float64, nscenarios, n)
    x = zeros(Float64, nscenarios, n)
    u = zeros(Float64, nscenarios, n)
    x_old = zeros(Float64, nscenarios, n)
    primres = dualres = Inf

    if hist!==nothing
        hist[:functionalvalue] = Float64[]
        hist[:computingtime] = Float64[]
        hist[:time] = Float64[]
        hist[:primres] = Float64[]
        hist[:dualres] = Float64[]
        hist[:residual] = Float64[]
        haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
        hist[:logstep] = printstep*nscenarios
    end

    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params

    # initialisation, for fair comparison with randomized methods.
    tinit = time()
    ph_initialization!(x, u, y, pb, μ, subpbparams, printlev)
    objval, dot_xu, primres, dualres = objective_value(pb, x), dot(pb, x, u), norm(pb, x-x_old), (1/μ) * norm(pb, y-x)

    it = 0
    computingtime = 0.0
    nscenariostreated = nscenarios

    printlev>0 && @printf " it   #scenario      primal res        dual res            dot(x,u)   objective\n"
    ph_print_log(pb, x, u, printlev, hist, it, nscenariostreated, primres, dualres, computingtime, tinit, callback)

    while (!ph_hasconverged(pb, x, primres, dualres, ϵ_abs, ϵ_rel) && it < maxiter
                                                    && time() - tinit < maxtime
                                                    && computingtime < maxcomputingtime)
        it += 1
        it_startcomputingtime = time()

        copy!(x_old, x)

        # Subproblem solves
        for id_scen in 1:nscenarios
            y[id_scen, :] = ph_subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, subpbparams)
        end
        nscenariostreated += n

        # projection on non anticipatory subspace
        nonanticipatory_projection!(x, pb, y)

        # multiplier update
        u += (1/μ) * (y-x)


        primres = norm(pb, x-x_old)
        dualres = (1/μ) * norm(pb, y-x)

        computingtime += time() - it_startcomputingtime

        # Print and logs
        if mod(it, printstep) == 0
            ph_print_log(pb, x, u, printlev, hist, it, nscenariostreated, primres, dualres, computingtime, tinit, callback)
        end
    end

    ## Final print
    if mod(it, printstep) != 0
        ph_print_log(pb, x, u, printlev, hist, it, nscenariostreated, primres, dualres, computingtime, tinit, callback)
    end

    printlev>0 && println("Computation time (s): ", computingtime)
    printlev>0 && println("Total time       (s): ", time() - tinit)

    return x
end
