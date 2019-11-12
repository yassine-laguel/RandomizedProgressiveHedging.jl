"""
    ph_subproblem_solve(pb, id_scen, u_scen, x_scen, μ, params)

Solve subproblem associated with scenario `id_scen`, primal point `x_scen`, dual point `u_scen`
and regularization `μ`.
"""
function ph_subproblem_solve(pb::Problem, id_scen::ScenarioId, u_scen, x_scen, μ, params::AbstractDict)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)

    # Augmented lagragian subproblem full objective
    obj += dot(u_scen, y) + (1/2*μ) * sum((y[i]-x_scen[i])^2 for i in 1:n)
    @objective(model, Min, obj)

    optimize!(model)
    if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
        @warn "Subproblem status: (scenario $id_scen) " primal_status(model) dual_status(model) termination_status(model)
    end

    return JuMP.value.(y)
end


"""
    x = solve_progressivehedging(pb::Problem)

Run the classical Progressive Hedging scheme on problem `pb`.

Stopping criterion is based on primal dual residual, maximum iterations or time
can also be set. Return a feasible point `x`.

## Keyword arguments:
- `ϵ_primal`: Tolerance on primal residual.
- `ϵ_dual`: Tolerance on dual residual.
- `μ`: Regularization parameter.
- `maxtime`: Limit time spent in computations.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
- `hist`: if not `nothing`, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
- `callback`: either nothing or a function `callback(pb, x, hist)::nothing` called at each log phase. `x` is the current feasible global iterate.
"""
function solve_progressivehedging(pb::Problem; ϵ_primal = 1e-4,
                                               ϵ_dual = 1e-4,
                                               μ = 3,
                                               maxtime = 3600,
                                               maxiter = 1e3,
                                               printlev = 1,
                                               printstep = 1,
                                               hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                               optimizer = Ipopt.Optimizer,
                                               optimizer_params = Dict{Symbol, Any}(:print_level=>0),
                                               callback = nothing,
                                               kwargs...)
    display_algopb_stats(pb, "Progressive Hedging", printlev, ϵ_primal=ϵ_primal, ϵ_dual=ϵ_dual, μ=μ, maxtime=maxtime, maxiter=maxiter)

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)

    y = zeros(Float64, nscenarios, n)
    x = zeros(Float64, nscenarios, n)
    u = zeros(Float64, nscenarios, n)
    primres = dualres = Inf

    (hist!==nothing) && (hist[:functionalvalue] = Float64[])
    (hist!==nothing) && (hist[:time] = Float64[])
    (hist!==nothing) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    (hist!==nothing) && (hist[:logstep] = printstep*nscenarios)

    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params

    it = 0
    tinit = time()
    printlev>0 && @printf " it   primal res        dual res            dot(x,u)   objective\n"
    while !(primres < ϵ_primal && dualres < ϵ_dual) && it < maxiter && time()-tinit < maxtime
        u_old = copy(u)

        # Subproblem solves
        for id_scen in 1:nscenarios
            y[id_scen, :] = ph_subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, subpbparams)
        end

        # projection on non anticipatory subspace
        nonanticipatory_projection!(x, pb, y)

        # multiplier update
        u += (1/μ) * (y-x)


        # invariants, indicators
        primres = norm(pb, x-y)
        dualres = (1/μ) * norm(pb, u - u_old)
        if mod(it, printstep) == 0
            objval = objective_value(pb, x)
            dot_xu = dot(pb, x, u)

            printlev>0 && @printf "%3i   %.10e  %.10e   % .3e % .16e\n" it primres dualres dot_xu objval

            (hist!==nothing) && push!(hist[:functionalvalue], objval)
            (hist!==nothing) && push!(hist[:time], time() - tinit)
            (hist!==nothing) && haskey(hist, :approxsol) && size(hist[:approxsol])==size(x) && push!(hist[:dist_opt], norm(hist[:approxsol] - x))

            if (callback!==nothing)
                callback(pb, x, hist)
            end
        end

        it += 1
    end

    ## Final print
    objval = objective_value(pb, x)
    dot_xu = dot(pb, x, u)

    printlev>0 && mod(it, printstep) != 1 && @printf "%3i   %.10e  %.10e   % .3e % .16e\n" it primres dualres dot_xu objval
    printlev>0 && println("Computation time: ", time() - tinit)

    return x
end
