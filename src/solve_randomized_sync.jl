"""
    randomizedsync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`xz_scen`)' where 'f_s' is the cost function associated with 
the scenario `id_scen`.
"""
function randomizedsync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params::AbstractDict)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)
    
    # Augmented lagragian subproblem full objective
    obj += (1/2*μ) * sum((y[i] - xz_scen[i])^2 for i in 1:n) ## TODO: replace with more explicit formula for efficiency
    @objective(model, Min, obj)
    
    optimize!(model)
    if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
        @warn "Subproblem of scenario $(id_scen) " primal_status(model) dual_status(model) termination_status(model)
    end

    return JuMP.value.(y)
end

"""
    randomizedsync_initialization!(z, pb, μ, subpbparams, printlev, it)

TODO
"""
function randomizedsync_initialization!(z, pb, μ, subpbparams, printlev, it)
    printlev>0 && print("Initialisation... ")
    
    xz_scen = zeros(get_scenariodim(pb))
    for id_scen in 1:pb.nscenarios
        z[id_scen, :] = randomizedsync_subpbsolve(pb, id_scen, xz_scen, μ, subpbparams)
        it += 1
    end
    nonanticipatory_projection!(z, pb, z)
    printlev>0 && println("done")
    return it
end

"""
    solve_randomized_sync(pb::Problem)

Run the Randomized Progressive Hedging scheme on problem `pb`.

Stopping criterion is maximum iterations or time. Return a feasible point `x`.

## Keyword arguments:
- `μ`: Regularization parameter.
- `qdistr`: if not nothing, specifies the probablility distribution for scenario sampling.
- `maxtime`: Limit time spent in computations.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
- `seed`: if not nothing, specifies the seed used for scenario sampling.
- `hist`: if not nothing, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
"""
function solve_randomized_sync(pb::Problem; μ = 3,
                                            qdistr = nothing,
                                            maxtime = 3600.0,
                                            maxiter = 1e4,
                                            printlev = 1,
                                            printstep = 1,
                                            seed = nothing,
                                            hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                            optimizer = Ipopt.Optimizer,
                                            optimizer_params = Dict{Symbol, Any}(:print_level=>0),
                                            kwargs...)
    printlev>0 && println("--------------------------------------------------------")
    printlev>0 && println("--- Randomized Progressive Hedging - synchronous")
    printlev>0 && println("--------------------------------------------------------")

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)
    
    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf
    
    rng = MersenneTwister(isnothing(seed) ? 1234 : seed)
    if isnothing(qdistr) || qdistr == :pdistr
        scen_sampling_distrib = Categorical(pb.probas)
    elseif qdistr == :unifdistr
        scen_sampling_distrib = Categorical(ones(nscenarios) / nscenarios)
    else
        @assert typeof(qdistr)<:Array
        @assert sum(qdistr) == 1
        @assert length(qdistr) == nscenarios
        scen_sampling_distrib = Categorical(qdistr)
    end

    !isnothing(hist) && (hist[:functionalvalue] = Float64[])
    !isnothing(hist) && (hist[:time] = Float64[])
    !isnothing(hist) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    !isnothing(hist) && (hist[:logstep] = printstep)

    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params

    it = 0
    tinit = time()
    printlev>0 && @printf "   it   global residual   objective\n"

    it = randomizedsync_initialization!(z, pb, μ, subpbparams, printlev, it)

    printlev>0 && @printf "%5i   %.10e % .16e\n" it 0.0 objective_value(pb, z)
    while it < maxiter && time()-tinit < maxtime
        id_scen = rand(rng, scen_sampling_distrib)

        ## Projection
        get_averagedtraj!(x, pb, z, id_scen) #TODO: rename with proj ?

        ## Subproblem solve
        y = randomizedsync_subpbsolve(pb, id_scen, 2*x-z[id_scen, :], μ, subpbparams)

        ## Global variable update
        z[id_scen, :] += (y - x)


        # invariants, indicators
        if mod(it, printstep) == 0
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            steplength = norm(x-y)
            
            printlev>0 && @printf "%5i   %.10e % .16e\n" it steplength objval

            !isnothing(hist) && push!(hist[:functionalvalue], objval)
            !isnothing(hist) && push!(hist[:time], time() - tinit)
            !isnothing(hist) && haskey(hist, :approxsol) && size(hist[:approxsol])==size(x_feas) && push!(hist[:dist_opt], norm(hist[:approxsol] - x_feas))
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
