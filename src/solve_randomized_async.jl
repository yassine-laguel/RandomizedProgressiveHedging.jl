struct AsyncSubproblemTask{T}
    taskid::Int
    scenario::T
    id_scenario::ScenarioId
    v_scen::Vector{Float64}
end

"""
    randasync_remote_func(inwork::RemoteChannel, outres::RemoteChannel, paramschan::RemoteChannel)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`v_scen`)' where 'f_s' is the cost function associated with
the scenario `id_scen`.
"""
function randasync_remote_func(inwork::RemoteChannel, outres::RemoteChannel, paramschan::RemoteChannel)
    try
        params = take!(paramschan)
        μ = params[:μ]
        build_fs = params[:build_fs]

        while true
            subpbtask::AsyncSubproblemTask = take!(inwork)

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

            put!(outres, (JuMP.value.(y), subpbtask.taskid))
        end
    catch e
        println("Worker error:")
        println(e)
        return
    end
end


"""
    randasync_init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario

Launch `asyncpar_remote_func` in every available worker with correct channels, and pass `subpbparams` to each worker.
"""
function randasync_init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario
    if workers() == Vector([1])
        @error "No workers available. Returning"
        return
    end

    nworkers = length(workers())

    work_channel = RemoteChannel(()->Channel{AsyncSubproblemTask{T}}(nworkers))
    results_channel = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Int}}(nworkers))
    params_channel = RemoteChannel(()->Channel{Dict{Symbol,Any}}(nworkers))

    remotecalls_futures = OrderedDict(worker_id => remotecall(randasync_remote_func,
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
    randasync_terminate_workers(pb, work_channel, remotecalls_futures)

Send a signal so that all workers return and wait that all workers do.
"""
function randasync_terminate_workers(pb, work_channel, remotecalls_futures)
    closing_task = AsyncSubproblemTask(
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
    randasync_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)

Compute a first global feasible point by solving each scenario independently and projecting
the global strategy obtained onto the non-anticipatory subspace. Independent solves are distributed
on available workers.
"""
function randasync_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)
    printlev>0 && print("Initialisation... ")
    xz_scen = zeros(get_scenariodim(pb))
    cur_scen = 1
    ## First, feed as many scenarios as there are workers
    for id_worker in 1:min(length(workers()), pb.nscenarios)
        put!(work_channel, AsyncSubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
        cur_scen += 1
    end

    it = 0
    while it < pb.nscenarios
        ## Taking current result
        y, id_scen = take!(results_channel)
        z[id_scen, :] = y

        ## Submitting new task if necessary
        if cur_scen <= pb.nscenarios
            put!(work_channel, AsyncSubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
            cur_scen += 1
        end
        it += 1
    end

    nonanticipatory_projection!(z, pb, z)

    printlev>0 && println("done")
    return
end


function randasync_init_hist!(hist, printstep)
    (hist!==nothing) && (hist[:functionalvalue] = Float64[])
    (hist!==nothing) && (hist[:computingtime] = Float64[])
    (hist!==nothing) && (hist[:time] = Float64[])
    (hist!==nothing) && (hist[:maxdelay] = Int32[])
    (hist!==nothing) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    (hist!==nothing) && (hist[:logstep] = printstep)
    return
end

function randasync_defaultsubpbparams(pb, optimizer, optimizer_params, μ)
    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params
    subpbparams[:μ] = μ
    subpbparams[:build_fs] = pb.build_subpb
    return subpbparams
end

function randasync_print_log(pb, z_feas, step, τ, delay, printlev, hist, it, computingtime, tinit, callback)
    objval = objective_value(pb, z_feas)
    steplength = norm(step)

    printlev>0 && @printf "%5i   %.10e   % .16e  %3i  %3i\n" it steplength objval τ delay

    (hist!==nothing) && push!(hist[:functionalvalue], objval)
    (hist!==nothing) && push!(hist[:computingtime], computingtime)
    (hist!==nothing) && push!(hist[:time], time() - tinit)
    (hist!==nothing) && haskey(hist, :approxsol) && size(hist[:approxsol])==size(z_feas) && push!(hist[:dist_opt], norm(hist[:approxsol] - z_feas))

    if (callback!==nothing)
        callback(pb, z_feas, hist)
    end
    return
end

"""
    solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario

Run the Randomized Progressive Hedging scheme on problem `pb`. All workers should be available.

Stopping criterion is maximum iterations or time. Return a feasible point `x`.

## Keyword arguments:
- `μ`: Regularization parameter.
- `c`: parameter for step length.
- `stepsize`: if `nothing` uses theoretical formula for stepsize, otherwise uses constant numerical value.
- `qdistr`: strategy of scenario sampling for subproblem selection,
    + `:pdistr`: samples according to objective scenario probability distribution,
    + `:unifdistr`: samples according to uniform distribution,
    + otherwise, an `Array` specifying directly the probablility distribution.
- `maxtime`: Limit on time spent in `solve_progressivehedging`.
- `maxcomputingtime`: Limit time spent in computations, excluding computation of initial feasible point and computations required by plottting / logging.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
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
function solve_randomized_async(pb::Problem{T}; μ::Float64 = 3.0,
                                                c = 0.9,
                                                stepsize::Union{Nothing, Float64} = nothing,
                                                qdistr = :pdistr,
                                                maxtime = 3600,
                                                maxcomputingtime = Inf,
                                                maxiter = 1e5,
                                                printlev = 1,
                                                printstep = 1,
                                                seed = nothing,
                                                hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                                optimizer = Ipopt.Optimizer,
                                                optimizer_params = Dict{String, Any}("print_level"=>0),
                                                callback=nothing,
                                                kwargs...) where T<:AbstractScenario
    display_algopb_stats(pb, "Randomized Progressive Hedging - asynchronous", printlev, μ=μ, c=c, stepsize=stepsize, qdistr=qdistr, maxtime=maxtime, maxcomputingtime=maxcomputingtime, maxiter=maxiter)

    # Variables
    nworkers = length(workers())
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)

    x = zeros(Float64, nworkers, n)
    z = zeros(Float64, nscenarios, n)
    z_feas = zeros(Float64, nscenarios, n)
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

    τ = ceil(Int, nworkers*1.05)                       # Upper bound on delay
    delay = 0

    randasync_init_hist!(hist, printstep)

    subpbparams = randasync_defaultsubpbparams(pb, optimizer, optimizer_params, μ)

    ## Workers Initialization
    printlev>0 && println("Available workers: ", nworkers)
    work_channel, results_channel, params_channel, remotecalls_futures = randasync_init_workers(pb, subpbparams)

    taskid_to_xcoord = Dict{Int, ScenarioId}()
    taskid_to_idscen = Dict{Int, ScenarioId}()
    taskid_to_launchit = Dict{Int, Int}()
    cur_taskid = 0

    # Algorithm initialisation
    tinit = time()
    randasync_initialization!(z, pb, μ, subpbparams, printlev, work_channel, results_channel)
    copy!(z_feas, z)

    it = 0
    computingtime = 0.0

    ## Feeding every worker with one task
    for (x_coord, w_id) in enumerate(workers())
        id_scen = rand(rng, scen_sampling_distrib)
        v = 2*x[x_coord, :] - z[id_scen, :]
        put!(work_channel, AsyncSubproblemTask(cur_taskid, pb.scenarios[id_scen], id_scen, v))

        taskid_to_xcoord[cur_taskid] = x_coord
        taskid_to_idscen[cur_taskid] = id_scen
        taskid_to_launchit[cur_taskid] = it
        cur_taskid += 1
    end

    objval = objective_value(pb, z)
    init_objval = objval

    printlev>0 && @printf "   it   residual            objective                 τ    delay\n"
    randasync_print_log(pb, z_feas, step, τ, delay, printlev, hist, it, computingtime, tinit, callback)

    while it < maxiter && time()-tinit < maxtime && computingtime < maxcomputingtime && objval < 3*init_objval
        it += 1
        it_startcomputingtime = time()

        ## Get a completed task
        y, taskid = take!(results_channel)
        x_coord = taskid_to_xcoord[taskid]; delete!(taskid_to_xcoord, taskid)
        id_scen = taskid_to_idscen[taskid]; delete!(taskid_to_idscen, taskid)
        delay = it-taskid_to_launchit[taskid]; delete!(taskid_to_launchit, taskid)

        ## z update
        τ = max(τ, delay)
        η = c*nscenarios*qmin / (2*τ*sqrt(qmin) + 1)

        if (stepsize === nothing)
            step = 2 * η / (nscenarios * pb.probas[id_scen]) * (y - x[x_coord, :])
        else
            step = stepsize * (y - x[x_coord, :])
        end
        z[id_scen, :] += step

        ## Draw new scenario, build v, send task
        id_scen = rand(scen_sampling_distrib)

        @views get_averagedtraj!(x[x_coord, :], pb, z, id_scen)
        v = 2*x[x_coord, :] - z[id_scen, :]

        put!(work_channel, AsyncSubproblemTask(cur_taskid, pb.scenarios[id_scen], id_scen, v))

        taskid_to_xcoord[cur_taskid] = x_coord
        taskid_to_idscen[cur_taskid] = id_scen
        taskid_to_launchit[cur_taskid] = it
        cur_taskid += 1


        computingtime += time() - it_startcomputingtime

        # Print and logs
        (hist!==nothing) && push!(hist[:maxdelay], delay)
        if mod(it, printstep) == 0
            nonanticipatory_projection!(z_feas, pb, z)
            randasync_print_log(pb, z_feas, step, τ, delay, printlev, hist, it, computingtime, tinit, callback)
        end

        it += 1
    end

    ## Get final solution
    nonanticipatory_projection!(z_feas, pb, z)
    randasync_print_log(pb, z_feas, step, τ, delay, printlev, hist, it, computingtime, tinit, callback)

    printlev>0 && println("Computation time (s): ", computingtime)
    printlev>0 && println("Total time       (s): ", time() - tinit)

    ## Terminate all workers
    randasync_terminate_workers(pb, work_channel, remotecalls_futures)

    return z_feas
end
