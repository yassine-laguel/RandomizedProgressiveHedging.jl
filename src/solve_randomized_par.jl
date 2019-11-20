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

            model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))

            # Get scenario objective function, build constraints in model
            y, obj, ctrref = build_fs(model, subpbtask.scenario, subpbtask.id_scenario)

            # Subproblem full objective
            obj += (1/(2*μ)) * sum((y[i] - subpbtask.v_scen[i])^2 for i in 1:length(y))

            @objective(model, Min, obj)

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
    randpar_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)

Compute a first global feasible point by solving each scenario independently and projecting 
the global strategy obtained onto the non-anticipatory subspace. Independent solves are distributed
on available workers.
"""
function randpar_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)
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
    return it
end


function randpar_init_hist!(hist, printstep)
    (hist!==nothing) && (hist[:functionalvalue] = Float64[])
    (hist!==nothing) && (hist[:time] = Float64[])
    (hist!==nothing) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    (hist!==nothing) && (hist[:logstep] = printstep)
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

"""
    solve_randomized_par(pb::Problem{T}) where T<:AbstractScenario

Run the Randomized Progressive Hedging scheme on problem `pb`. All workers should be available.

Stopping criterion is maximum iterations or time. Return a feasible point `x`.

## Keyword arguments:
- `μ`: Regularization parameter.
- `c`: parameter for step length.
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
    + `:number_waitingworkers`: array of number of wainting workers, indexed by iteration,
    + `:maxdelay`: array of maximal delay among done workers, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
- `callback`: either `nothing` or a function `callback(pb, x, hist)::nothing` called at each 
log phase. `x` is the current feasible global iterate.
"""
function solve_randomized_par(pb::Problem{T}; μ::Float64 = 3.0,
                                                c = 0.9,
                                                qdistr = :pdistr,
                                                maxtime = 3600,
                                                maxiter = 1e5,
                                                printlev = 1,
                                                printstep = 1,
                                                seed = nothing,
                                                hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                                optimizer = Ipopt.Optimizer,
                                                optimizer_params = Dict{Symbol, Any}(:print_level=>0),
                                                callback=nothing,
                                                kwargs...) where T<:AbstractScenario
    display_algopb_stats(pb, "Randomized Progressive Hedging - parallel", printlev, μ=μ, c=c, qdistr=qdistr, maxtime=maxtime, maxiter=maxiter)

    # Variables
    nworkers = length(workers())
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)

    x = zeros(Float64, nworkers, n)
    z = zeros(Float64, nscenarios, n)
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
    work_channel, results_channel, params_channel, remotecalls_futures = randpar_init_workers(pb, subpbparams)

    it = 0
    tinit = time()
    printlev>0 && @printf "   it   residual            objective                \n"

    it = randpar_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)
    x = copy(z)

    printlev>0 && @printf "%5i   %.10e   % .16e    \n" it 0.0 objective_value(pb, z)


    minworkscen = min(nworkers, nscenarios)
    minworkscen < nworkers && @warn "solve_randomized_par(): less scenarios than workers ($nscenarios < $nworkers). Only $nscenarios workers will be put to use."


    while it < maxiter && time()-tinit < maxtime

        z_prev = copy(z)

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

        x = nonanticipatory_projection(pb, z)

        if mod(it, printstep) == 0
            objval = objective_value(pb, x)
            steplength = norm(z-z_prev)

            printlev>0 && @printf "%5i   %.10e   % .16e \n" it steplength objval

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
    steplength = norm(z-z_prev)

    printlev>0 && mod(it, printstep) != 1 && @printf "%5i   %.10e   % .16e\n" it steplength objval
    printlev>0 && println("Computation time: ", time() - tinit)

    ## Terminate all workers
    randpar_terminate_workers(pb, work_channel, remotecalls_futures)

    return x
end
