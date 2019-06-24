struct SubproblemTask{T}
    scenario::T
    id_scenario::ScenarioId
    build_subpb::Function
    μ::Float64
    v_scen::Vector{Float64}
end

"""
do_remote_work(inwork::RemoteChannel, outres::RemoteChannel)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`v_scen`)' where 'f_s' is the cost function associated with 
the scenario `id_scen`.
"""
function do_remote_work(inwork::RemoteChannel, outres::RemoteChannel, paramschan::RemoteChannel)
    params = take!(paramschan)
    while true
        try
            t0 = time()
            subpbtask::SubproblemTask = take!(inwork)

            if subpbtask.id_scenario == -1  # Work finished
                return
            end

            # do work
            model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))
            
            # Get scenario objective function, build constraints in model
            y, obj, ctrref = subpbtask.build_subpb(model, subpbtask.scenario, subpbtask.id_scenario)
            
            obj += (1/2*subpbtask.μ) * sum((y[i] - subpbtask.v_scen[i])^2 for i in 1:length(y))
            
            @objective(model, Min, obj)

            optimize!(model)

            put!(outres, JuMP.value.(y))        
        catch e
            println("Worker error:")
            println(e)
        end
    end
end



function init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario
    if workers() == Vector([1])
        @error "No workers available. Returning"
        return
    end

    work_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{SubproblemTask{T}}(3), worker_id) for worker_id in workers())
    results_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Vector{Float64}}(3), worker_id) for worker_id in workers())
    params_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Dict{Symbol,Any}}(3), worker_id) for worker_id in workers())
    
    remotecalls_futures = OrderedDict(worker_id => remotecall(do_remote_work, worker_id, work_channels[worker_id], results_channels[worker_id], params_channels[worker_id]) for worker_id in workers())
    for w_id in workers()
        put!(params_channels[w_id], subpbparams)
    end

    return work_channels, results_channels, params_channels, remotecalls_futures
end

function terminate_workers(pb, work_channels, remotecalls_futures)
    closing_task = SubproblemTask(
        pb.scenarios[1],
        -1,                 ## Stop signal
        isnothing,
        0.0,
        [0.0],
    )
    for (w_id, worker_wkchan) in work_channels
        put!(worker_wkchan, closing_task)
    end

    for (w_id, worker) in remotecalls_futures
        wait(worker)
    end
    return
end

function get_oldest_readyworker(results_channels, it, worker_to_launchit, nwaitingworkers, maxdelay)
    ## Wait for a worker to complete its job
    ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
    while length(ready_workers) == 0
        ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
    end
    
    ## Get oldest ready worker id
    worker_to_delay = SortedDict(w_id=>it-worker_to_launchit[w_id] for w_id in ready_workers)

    return argmax(worker_to_delay), max(nwaitingworkers, length(ready_workers)), max(maxdelay, maximum(values(worker_to_delay)))
end

function giveworkertask!(work_channels, worker_to_launchit, worker_to_scen, it, cur_worker, pb, id_scen, μ, v)
    task = SubproblemTask(
        pb.scenarios[id_scen],
        id_scen,
        pb.build_subpb,
        μ,
        v,
    )
    
    put!(work_channels[cur_worker], task)
    worker_to_launchit[cur_worker] = it
    worker_to_scen[cur_worker] = id_scen
    return
end


"""
    solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario

Run the Randomized Progressive Hedging scheme on problem `pb`. All workers should be available.

Stopping criterion is maximum iterations or time. Return a feasible point `x`.

## Keyword arguments:
- `μ`: Regularization parameter.
- `c`: parameter for step length.
- `qdistr`: if not nothing, specifies the probablility distribution for scenario sampling.
- `maxtime`: Limit time spent in computations.
- `maxiter`: Maximum iterations.
- `printlev`: if 0, mutes output.
- `printstep`: number of iterations skipped between each print and logging.
- `seed`: if not nothing, specifies the seed used for scenario sampling.
- `hist`: if not nothing, will record:
    + `:functionalvalue`: array of functional value indexed by iteration,
    + `:time`: array of time at the end of each iteration, indexed by iteration,
    + `:number_waitingworkers`: array of number of wainting workers, indexed by iteration,
    + `:maxdelay`: array of maximal delay among done workers, indexed by iteration,
    + `:dist_opt`: if dict has entry `:approxsol`, array of distance between iterate and `hist[:approxsol]`, indexed by iteration.
    + `:logstep`: number of iteration between each log.
- `optimizer`: an optimizer for subproblem solve.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
"""
function solve_randomized_async(pb::Problem{T}; μ::Float64 = 3.0,
                                                c = 0.9,
                                                qdistr = nothing,
                                                maxtime = 3600,
                                                maxiter = 1e5,
                                                printlev = 1,
                                                printstep = 1,
                                                seed = nothing,
                                                hist::Union{OrderedDict{Symbol, Any}, Nothing}=nothing,
                                                optimizer = Ipopt.Optimizer,
                                                optimizer_params = Dict{Symbol, Any}(:print_level=>0),
                                                kwargs...) where T<:AbstractScenario
    printlev>0 && println("--------------------------------------------------------")
    printlev>0 && println("--- Randomized Progressive Hedging - asynchronous")
    printlev>0 && println("--------------------------------------------------------")
    
    # Variables
    nworkers = length(workers())
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))
    
    x = zeros(Float64, nworkers, n)
    step = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf
    
    rng = MersenneTwister(isnothing(seed) ? 1234 : seed)
    scen_sampling_distrib = Categorical(isnothing(qdistr) ? pb.probas : qdistr)
    qmin = minimum(scen_sampling_distrib.p)      
    τ = ceil(Int, nworkers*1.05)                       # Upper bound on delay
    
    nwaitingworkers = maxdelay = 0
    !isnothing(hist) && (hist[:functionalvalue] = Float64[])
    !isnothing(hist) && (hist[:time] = Float64[])
    !isnothing(hist) && (hist[:number_waitingworkers] = Int32[])
    !isnothing(hist) && (hist[:maxdelay] = Int32[])
    !isnothing(hist) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    !isnothing(hist) && (hist[:logstep] = printstep)

    subpbparams = OrderedDict{Symbol, Any}()
    subpbparams[:optimizer] = optimizer
    subpbparams[:optimizer_params] = optimizer_params

    ## Workers Initialization
    printlev>0 && println("Available workers: ", nworkers)

    printlev>0 && printstyled("Launching remotecalls... ", color=:red)
    work_channels, results_channels, params_channels, remotecalls_futures = init_workers(pb, subpbparams)
    printlev>0 && printstyled("Done.\n", color=:red)
    worker_to_scen = OrderedDict{Int, ScenarioId}()
    worker_to_launchit = OrderedDict{Int, ScenarioId}()
    worker_to_wid = SortedDict{Int, Int}(worker=>wid for (wid, worker) in enumerate(workers())) # workers() -> 1:nworkers map for x matrix

    ## Feeding every worker with one task
    for w_id in workers()
        id_scen = rand(rng, scen_sampling_distrib)
        giveworkertask!(work_channels, worker_to_launchit, worker_to_scen, 0, w_id, pb, id_scen, μ, zeros(n))
    end

    it = 0
    tinit = time()
    printlev>0 && @printf "   it   residual            objective                 τ    maxdelay  #waitingworkers\n"
    while it < maxiter && time()-tinit < maxtime
        ## Wait for a worker to complete job, take max delay one if several
        cur_worker, nwaitingworkers, maxdelay = get_oldest_readyworker(results_channels, it, worker_to_launchit, nwaitingworkers, maxdelay)

        y = take!(results_channels[cur_worker])
        id_scen = worker_to_scen[cur_worker]

        ## z update
        τ = max(τ, it-minimum(values(worker_to_launchit)))
        η = c*nscenarios*qmin / (2*τ*sqrt(qmin) + 1)
    
        step = 2 * η / (nscenarios * pb.probas[id_scen]) * (y - x[worker_to_wid[cur_worker], :])
        z[id_scen, :] += step
        
        ## Draw new scenario for worker, build v and task
        id_scen = rand(scen_sampling_distrib)
        
        @views get_averagedtraj!(x[worker_to_wid[cur_worker], :], pb, z, id_scen)
        v = 2*x[worker_to_wid[cur_worker], :] - z[id_scen, :]

        giveworkertask!(work_channels, worker_to_launchit, worker_to_scen, it, cur_worker, pb, id_scen, μ, v)

        # invariants, indicators, prints
        !isnothing(hist) && push!(hist[:number_waitingworkers], nwaitingworkers)
        !isnothing(hist) && push!(hist[:maxdelay], maxdelay)

        if mod(it, printstep) == 0
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            steplength = norm(step)
            
            printlev>0 && @printf "%5i   %.10e   % .16e  %3i  %3i       %3i\n" it steplength objval τ maxdelay nwaitingworkers

            !isnothing(hist) && push!(hist[:functionalvalue], objval)
            !isnothing(hist) && push!(hist[:time], time() - tinit)
            !isnothing(hist) && haskey(hist, :approxsol) && push!(hist[:dist_opt], norm(hist[:approxsol] - x_feas))
            
            nwaitingworkers = maxdelay = 0
        end
        
        it += 1
    end
    
    ## Get final solution
    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    
    printlev>0 && mod(it, printstep) == 1 && @printf "%5i   %.10e   % .16e  %3i  %3i       %3i\n" it steplength objval τ maxdelay nwaitingworkers
    printlev>0 && println("Computation time: ", time() - tinit)
    
    ## Terminate all workers
    printlev>0 && printstyled("Terminating nodes...\n", color=:red)
    terminate_workers(pb, work_channels, remotecalls_futures)

    return x_feas
end
