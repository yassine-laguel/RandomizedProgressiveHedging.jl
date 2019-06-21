"""
randomizedsync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`xz_scen`)' where 'f_s' is the cost function associated with 
the scenario `id_scen`.
"""
function randomizedsync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params::AbstractDict)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=get(params, :print_level, 0)))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)
    
    # Augmented lagragian subproblem full objective
    obj += (1/2*μ) * sum((y[i] - xz_scen[i])^2 for i in 1:n) ## TODO: replace with more explicit formula for efficiency
    @objective(model, Min, obj)
    
    optimize!(model)
    
    return JuMP.value.(y)
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
"""
function solve_randomized_sync(pb::Problem; μ = 3,
                                            qdistr = nothing,
                                            maxtime = 3600,
                                            maxiter = 1e4,
                                            printlev = 1,
                                            printstep = 1,
                                            seed = nothing,
                                            hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                            kwargs...)
    printlev>0 && println("--------------------------------------------------------")
    printlev>0 && println("--- Randomized Progressive Hedging - synchronous")
    printlev>0 && println("--------------------------------------------------------")

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))
    
    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf
    
    rng = MersenneTwister(isnothing(seed) ? 1234 : seed)
    scen_sampling_distrib = Categorical(isnothing(qdistr) ? pb.probas : qdistr)
    
    !isnothing(hist) && (hist[:functionalvalue] = Float64[])
    !isnothing(hist) && (hist[:time] = Float64[])
    !isnothing(hist) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    !isnothing(hist) && (hist[:logstep] = printstep)

    ## Subpb
    subpbparams = OrderedDict()

    it = 0
    tinit = time()
    printlev>0 && @printf "   it   global residual   objective\n"
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
            !isnothing(hist) && haskey(hist, :approxsol) && push!(hist[:dist_opt], norm(hist[:approxsol] - x_feas))
        end
        
        it += 1
    end

    ## Final print
    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    steplength = norm(x-y)
    
    printlev>0 && mod(it, printstep) == 1 && @printf "%5i   %.10e % .16e\n" it steplength objval
    printlev>0 && println("Computation time: ", time() - tinit)

    return x_feas
end
