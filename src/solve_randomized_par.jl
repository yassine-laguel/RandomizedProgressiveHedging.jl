struct ParSubproblemTask{T}
    taskid::Int
    scenario::T
    id_scenario::ScenarioId
    v_scen::Vector{Float64}
end

"""
    randpar_remote_func(inwork::RemoteChannel, outres::RemoteChannel)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`v_scen`)' where 'f_s' is the cost function associated with
the scenario `id_scen`.
"""
function randpar_remote_func(inwork::RemoteChannel, outres::RemoteChannel, paramschan::RemoteChannel)
    try
        params = take!(paramschan)
        μ = params[:μ]
        build_fs = params[:build_fs]

        while true
            subpbtask::ParSubproblemTask = take!(inwork)

            if subpbtask.taskid == -1  # Work finished
                return
            end

            model = Model(optimizer_with_attributes(params[:optimizer], params[:optimizer_params]...))

            # Get scenario objective function, build constraints in model
            y, obj, ctrref = build_fs(model, subpbtask.scenario, subpbtask.id_scenario)

            # Subproblem full objective
            obj += (1/(2*μ)) * sum((y[i] - subpbtask.v_scen[i])^2 for i in 1:length(y))

            @objective(model, MOI.MIN_SENSE, obj)

            optimize!(model)
            if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
                @warn "Subproblem of scenario $(subpbtask.id_scenario) " primal_status(model) dual_status(model) termination_status(model)
            end

            put!(outres, (JuMP.value.(y), subpbtask.id_scenario))
        end
    catch e
        println("Worker error:")
        println(e)
        return
    end
end


"""
    randpar_init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario

Launch `randpar_remote_func` in every available worker with correct channels, and pass `subpbparams` to each worker.
"""
function randpar_init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario
    if workers() == Vector([1])
        @error "No workers available. Returning"
        return
    end

    nworkers = length(workers())

    work_channel = RemoteChannel(()->Channel{ParSubproblemTask{T}}(nworkers))
    results_channel = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Int}}(nworkers))
    params_channel = RemoteChannel(()->Channel{Dict{Symbol,Any}}(nworkers))

    remotecalls_futures = OrderedDict(worker_id => remotecall(randpar_remote_func,
                                                              worker_id,
                                                              work_channel,
                                                              results_channel,
                                                              params_channel) for worker_id in workers())
    for w_id in workers()
        put!(params_channel, subpbparams)
    end

    return work_channel, results_channel, params_channel, remotecalls_futures
end

"""
    randpar_terminate_workers(pb, work_channel, remotecalls_futures)

Send a signal so that all workers return and wait that all workers do.
"""
function randpar_terminate_workers(pb, work_channel, remotecalls_futures)
    closing_task = ParSubproblemTask(
        -1,                 ## Stop signal
        pb.scenarios[1],
        0,
        [0.0]
    )
    for w_id in workers()
        put!(work_channel, closing_task)
    end

    for (w_id, worker) in remotecalls_futures
        wait(worker)
    end
    return
end

"""
    randpar_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)

Compute a first global feasible point by solving each scenario independently and projecting
the global strategy obtained onto the non-anticipatory subspace. Independent solves are distributed
on available workers.
"""
function randpar_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)
    printlev>0 && print("Initialisation... ")
    xz_scen = zeros(get_scenariodim(pb))
    cur_scen = 1
    ## First, feed as many scenarios as there are workers
    for id_worker in 1:min(length(workers()), pb.nscenarios)
        put!(work_channel, ParSubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
        cur_scen += 1
    end

    it = 0
    while it < pb.nscenarios
        ## Taking current result
        y, id_scen = take!(results_channel)
        z[id_scen, :] = y

        ## Submitting new task if necessary
        if cur_scen <= pb.nscenarios
            put!(work_channel, ParSubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
            cur_scen += 1
        end
        it += 1
    end

    nonanticipatory_projection!(z, pb, z)

    printlev>0 && println("done")
    return
end


function randpar_init_hist!(hist, printstep)
    if hist!==nothing
        hist[:functionalvalue] = Float64[]
        hist[:residual] = Float64[]
        hist[:computingtime] = Float64[]
        hist[:time] = Float64[]
        hist[:logstep] = printstep
        haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    end
    return
end

function randpar_defaultsubpbparams(pb, optimizer, optimizer_params, μ)
    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params
    subpbparams[:μ] = μ
    subpbparams[:build_fs] = pb.build_subpb
    return subpbparams
end

function randpar_print_log(pb, z, z_prev, x, printlev, residual, hist, it, nscenariostreated, computingtime, tinit, callback)
    objval = objective_value(pb, x)
    steplength = norm(z-z_prev)

    printlev>0 && @printf "%5i   %.2e       %.10e   %.10e     % .16e\n" it nscenariostreated steplength residual objval

    if hist!==nothing
        push!(hist[:functionalvalue], objval)
        push!(hist[:residual], residual)
        push!(hist[:computingtime], computingtime)
        push!(hist[:time], time() - tinit)
        haskey(hist, :approxsol) && size(hist[:approxsol])==size(x) && push!(hist[:dist_opt], norm(hist[:approxsol] - x))
    end

    if (callback!==nothing)
        callback(pb, x, hist)
    end
    return
end

"""
    solve_randomized_par(pb::Problem{T}) where T<:AbstractScenario

Run the Randomized Progressive Hedging scheme on problem `pb`. All workers should be available.

Stopping criterion is maximum iterations, time, or the residual, defined as ``|z_{k-d}-z_k\``
where ``d`` is the interval between checks, `residualupdate_interval`, equal to the number
of scenarios by default. Return a feasible point `x`.

## Keyword arguments:
- `ϵ_abs`: Absolute tolerance on residual.
- `ϵ_rel`: Relative tolerance on residual.
- `residualupdate_interval`: number of scenario to be treated before computing and checking
residual, default value is the number of scenarios in the problem.
- `μ`: Regularization parameter.
- `c`: parameter for step length.
- `qdistr`: strategy of scenario sampling for subproblem selection,
    + `:pdistr`: samples according to objective scenario probability distribution,
    + `:unifdistr`: samples according to uniform distribution,
    + otherwise, an `Array` specifying directly the probablility distribution.
- `maxtime`: Limit on time spent in `solve_progressivehedging`.
- `maxcomputingtime`: Limit time spent in computations, excluding computation of initial feasible point and computations required by plottting / logging.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging, default value is
the number of scenarios in problem.
- `seed`: if not `nothing`, specifies the seed used for scenario sampling.
- `hist`: if not `nothing`, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:number_waitingworkers`: array of number of wainting workers, indexed by iteration,
    + `:maxdelay`: array of maximal delay among done workers, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{String, Any}` storing parameters for the optimizer.
- `callback`: either `nothing` or a function `callback(pb, x, hist)::nothing` called at each
log phase. `x` is the current feasible global iterate.
"""
function solve_randomized_par(pb::Problem{T}; ϵ_abs = 1e-8,
                                              ϵ_rel = 1e-4,
                                              residualupdate_interval = nothing,
                                              μ::Float64 = 3.0,
                                              c = 0.9,
                                              qdistr = :pdistr,
                                              maxtime = 3600,
                                              maxcomputingtime = Inf,
                                              maxiter = 1e5,
                                              printlev = 1,
                                              printstep = nothing,
                                              seed = nothing,
                                              hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                              optimizer = Ipopt.Optimizer,
                                              optimizer_params = Dict{String, Any}("print_level"=>0),
                                              callback=nothing,
                                              kwargs...) where T<:AbstractScenario
    display_algopb_stats(pb, "Randomized Progressive Hedging - parallel", printlev, ϵ_abs=ϵ_abs, ϵ_rel=ϵ_rel, μ=μ, c=c, qdistr=qdistr, maxtime=maxtime, maxcomputingtime=maxcomputingtime, maxiter=maxiter)

    # Variables
    nworkers = length(workers())
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)
    isnothing(residualupdate_interval) && (residualupdate_interval = nscenarios)
    isnothing(printstep) && (printstep = ceil(nscenarios / min(nworkers, nscenarios)))

    x = zeros(Float64, nscenarios, n)
    z = zeros(Float64, nscenarios, n)
    z_old = zeros(Float64, nscenarios, n)
    z_prev = zeros(Float64, nscenarios, n)
    step = zeros(Float64, n)
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
    qmin = minimum(scen_sampling_distrib.p)


    randpar_init_hist!(hist, printstep*min(nworkers,pb.nscenarios))

    subpbparams = randpar_defaultsubpbparams(pb, optimizer, optimizer_params, μ)

    ## Workers Initialization
    printlev>0 && println("Available workers: ", nworkers)
    minworkscen = min(nworkers, nscenarios)
    minworkscen < nworkers && @warn "solve_randomized_par(): less scenarios than workers ($nscenarios < $nworkers). Only $nscenarios workers will be put to use."

    work_channel, results_channel, params_channel, remotecalls_futures = randpar_init_workers(pb, subpbparams)

    # Algorithm initialisation
    tinit = time()
    randpar_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)
    copy!(x, z)

    it = 0
    computingtime = 0.0
    nscenariostreated = nscenarios
    lastresidualupdate = 0
    residual = rand_getresidual(pb, z, z_old)

    printlev>0 && @printf "   it   #scenarios     iteration step     residual              objective\n"
    randpar_print_log(pb, z, z_prev, x, printlev, residual, hist, it, nscenariostreated, computingtime, tinit, callback)

    while !rand_hasconverged(pb, z, residual, ϵ_abs, ϵ_rel) &&
            it < maxiter &&
            time()-tinit < maxtime &&
            computingtime < maxcomputingtime
        it += 1
        nscenariostreated += minworkscen
        it_startcomputingtime = time()

        copy!(z_prev, z)

        ## Draw a scenario, build v, send task
        tab_id_scen = zeros(Int64, minworkscen)
        for i in 1:minworkscen
            id_scen = rand(rng, scen_sampling_distrib)
            while id_scen in tab_id_scen[1:i-1]
                id_scen = rand(rng, scen_sampling_distrib)
            end
            tab_id_scen[i] = id_scen
        end

        ## For all available workers
        for id_scen in tab_id_scen
            v = 2*x[id_scen, :] - z[id_scen, :]
            put!(work_channel, ParSubproblemTask(id_scen, pb.scenarios[id_scen], id_scen, v))
        end

        for resp_worker in tab_id_scen
            ## Taking current result
            y, id_scen = take!(results_channel)

            z[id_scen, :] += (y - x[id_scen, :])
        end

        nonanticipatory_projection!(x, pb, z)


        computingtime += time() - it_startcomputingtime


        # update residual if enough scenario have been treated since last residual eval
        if nscenariostreated > lastresidualupdate + residualupdate_interval
            residual = rand_getresidual(pb, z, z_old)
            copy!(z_old, z)
            lastresidualupdate = nscenariostreated
        end

        # Print and logs
        if mod(it, printstep) == 0
            randpar_print_log(pb, z, z_prev, x, printlev, residual, hist, it, nscenariostreated, computingtime, tinit, callback)
        end
    end

    ## Final print
    if mod(it, printstep) != 0
        randpar_print_log(pb, z, z_prev, x, printlev, residual, hist, it, nscenariostreated, computingtime, tinit, callback)
    end

    printlev>0 && println("Computation time (s): ", computingtime)
    printlev>0 && println("Total time       (s): ", time() - tinit)

    ## Terminate all workers
    randpar_terminate_workers(pb, work_channel, remotecalls_futures)

    return x
end
