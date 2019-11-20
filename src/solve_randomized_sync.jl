"""
    randsync_subproblem_solve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`xz_scen`)' where 'f_s' is the cost function associated with
the scenario `id_scen`.
"""
function randsync_subproblem_solve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params::AbstractDict)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)

    # Subproblem full objective
    obj += (1/(2*μ)) * sum((y[i] - xz_scen[i])^2 for i in 1:n)
    @objective(model, Min, obj)

    optimize!(model)
    if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
        @warn "Subproblem of scenario $(id_scen) " primal_status(model) dual_status(model) termination_status(model)
    end

    return JuMP.value.(y)
end

"""
    randsync_initialization!(z, pb, μ, subpbparams, printlev, it)

Compute a first global feasible point by solving each scenario independently and projecting 
the global strategy obtained onto the non-anticipatory subspace.
"""
function randsync_initialization!(z, pb, μ, subpbparams, printlev, it)
    printlev>0 && print("Initialisation... ")

    xz_scen = zeros(get_scenariodim(pb))
    for id_scen in 1:pb.nscenarios
        z[id_scen, :] = randsync_subproblem_solve(pb, id_scen, xz_scen, μ, subpbparams)
        it += 1
    end
    nonanticipatory_projection!(z, pb, z)
    printlev>0 && println("done")
    return it
end

"""
    x = solve_randomized_sync(pb::Problem)

Run the Randomized Progressive Hedging scheme on problem `pb`.

Stopping criterion is maximum iterations or time. Return a feasible point `x`.

## Keyword arguments:
- `μ`: Regularization parameter.
- `qdistr`: strategy of scenario sampling for subproblem selection,
    + `:pdistr`: samples according to objective scenario probability distribution,
    + `:unifdistr`: samples according to uniform distribution,
    + otherwise, an `Array` specifying directly the probablility distribution.
- `maxtime`: Limit time spent in computations.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
- `seed`: if not `nothing`, specifies the seed used for scenario sampling.
- `hist`: if not `nothing`, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and 
    `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
- `callback`: either `nothing` or a function `callback(pb, x, hist)::nothing` called at each 
log phase. `x` is the current feasible global iterate.
"""
function solve_randomized_sync(pb::Problem; μ::Float64 = 3.0,
                                            qdistr = :pdistr,
                                            maxtime = 3600.0,
                                            maxiter = 1e4,
                                            printlev = 1,
                                            printstep = 1,
                                            seed = nothing,
                                            hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                            optimizer = Ipopt.Optimizer,
                                            optimizer_params = Dict{Symbol, Any}(:print_level=>0),
                                            callback=nothing,
                                            kwargs...)
    display_algopb_stats(pb, "Randomized Progressive Hedging - synchronous", printlev, μ=μ, qdistr=qdistr, maxtime=maxtime, maxiter=maxiter)

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)

    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf

    ## Random scenario sampling
    rng = MersenneTwister(seed===nothing ? 1234 : seed)
    if (qdistr===nothing) || qdistr == :pdistr
        scen_sampling_distrib = Categorical(pb.probas)
    elseif qdistr == :unifdistr
        scen_sampling_distrib = Categorical(ones(nscenarios) / nscenarios)
    else
        @assert typeof(qdistr)<:Array
        @assert sum(qdistr) == 1
        @assert length(qdistr) == nscenarios
        scen_sampling_distrib = Categorical(qdistr)
    end

    (hist!==nothing) && (hist[:functionalvalue] = Float64[])
    (hist!==nothing) && (hist[:time] = Float64[])
    (hist!==nothing) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    (hist!==nothing) && (hist[:logstep] = printstep)

    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params

    it = 0
    tinit = time()
    printlev>0 && @printf "   it   global residual   objective\n"

    it = randsync_initialization!(z, pb, μ, subpbparams, printlev, it)

    printlev>0 && @printf "%5i   %.10e % .16e\n" it 0.0 objective_value(pb, z)
    while it < maxiter && time()-tinit < maxtime
        id_scen = rand(rng, scen_sampling_distrib)

        ## Projection
        get_averagedtraj!(x, pb, z, id_scen)

        ## Subproblem solve
        y = randsync_subproblem_solve(pb, id_scen, 2*x-z[id_scen, :], μ, subpbparams)

        ## Global variable update
        z[id_scen, :] += (y - x)


        # invariants, indicators
        if mod(it, printstep) == 0
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            steplength = norm(x-y)

            printlev>0 && @printf "%5i   %.10e % .16e\n" it steplength objval

            (hist!==nothing) && push!(hist[:functionalvalue], objval)
            (hist!==nothing) && push!(hist[:time], time() - tinit)
            (hist!==nothing) && haskey(hist, :approxsol) && size(hist[:approxsol])==size(x_feas) && push!(hist[:dist_opt], norm(hist[:approxsol] - x_feas))

            if (callback!==nothing)
                callback(pb, x_feas, hist)
            end

        end

        it += 1
    end

    ## Final print
    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    steplength = norm(x-y)

    printlev>0 && mod(it, printstep) != 1 && @printf "%5i   %.10e % .16e\n" it steplength objval
    printlev>0 && println("Computation time: ", time() - tinit)

    return x_feas
end
