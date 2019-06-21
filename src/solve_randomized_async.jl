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
function do_remote_work(inwork::RemoteChannel, outres::RemoteChannel)
    while true
        try
            t0 = time()
            subpbtask::SubproblemTask = take!(inwork)

            if subpbtask.id_scenario == -1  # Work finished
                return
            end

            # do work
            model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
            
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



function init_workers(pb::Problem{T}) where T<:AbstractScenario
    work_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{SubproblemTask{T}}(3), worker_id) for worker_id in workers())
    results_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Vector{Float64}}(3), worker_id) for worker_id in workers())
    
    printstyled("Launching remotecalls... ", color=:red)
    remotecalls_futures = OrderedDict(worker_id => remotecall(do_remote_work, worker_id, work_channels[worker_id], results_channels[worker_id]) for worker_id in workers())
    printstyled("Done.\n", color=:red)
    return work_channels, results_channels, remotecalls_futures
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

"""
    solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario

Run the Randomized Progressive Hedging scheme on problem `pb`. All workers should be available.
"""
function solve_randomized_async(pb::Problem{T}, kwargs...) where T<:AbstractScenario
    println("--------------------------------------------------------")
    println("--- Randomized Progressive Hedging - asynchronous")
    println("--------------------------------------------------------")
    println()
    println("Available workers: ", workers())
    nworkers = length(workers())

    if workers() == Vector([1])
        @error "No workers available. Returning"
        return zeros(1, 1).-1
    end
    
    # parameters
    μ = 3.0
    params = Dict(
        :print_step=>10,
        :max_iter=>800,
        :step_tol=>1e-4,
    )
    
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))
    
    qproba = ones(nscenarios) / nscenarios
    scen_sampling_distrib = Categorical(qproba)

    c = 0.9
    qmin = minimum(qproba)      
    τ = ceil(Int, nworkers*1.1)                       # Upper bound on delay
    η = c*nscenarios*qmin / (2*τ*sqrt(qmin) + 1)

    # Optim. variables
    x = zeros(Float64, nworkers, n)
    step = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf

    t_init = time()
    

    # Variables Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    ## Workers Initialization
    work_channels, results_channels, remotecalls_futures = init_workers(pb)

    worker_to_scen = OrderedDict{Int, ScenarioId}()
    worker_to_launchit = OrderedDict{Int, ScenarioId}()
    worker_to_wid = SortedDict{Int, Int}(worker=>wid for (wid, worker) in enumerate(workers())) # workers() -> 1:nworkers map for x matrix

    ## Feeding every worker with one task
    for w_id in workers()
        ## Draw random scenario, build subproblem task
        id_scen = rand(scen_sampling_distrib)

        task = SubproblemTask(
            pb.scenarios[id_scen],
            id_scen,
            pb.build_subpb,
            μ,
            zeros(n),
        )

        worker_to_scen[w_id] = id_scen
        worker_to_launchit[w_id] = 0

        put!(work_channels[w_id], task)
    end


    it = 0
    @printf "\n it    global residual   objective\n"
    while !(steplength < params[:step_tol]) && it < params[:max_iter]

        ## Wait for a worker to complete its job
        ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
        while length(ready_workers) == 0
            ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
        end
        
        ## Get oldest worker id
        worker_to_delay = SortedDict(w_id=>it-worker_to_launchit[w_id] for w_id in ready_workers)
        cur_worker = argmax(worker_to_delay)
        
        ## Collect result
        y = take!(results_channels[cur_worker])
        id_scen = worker_to_scen[cur_worker]

        ## z update
        step = 2 * η / (nscenarios * pb.probas[id_scen]) * (y - x[worker_to_wid[cur_worker], :])
        z[id_scen, :] += step
        
        ## Draw new scenario for worker, build v and task
        id_scen = rand(scen_sampling_distrib)
        
        @views get_averagedtraj!(x[worker_to_wid[cur_worker], :], pb, z, id_scen)
        v = 2*x[worker_to_wid[cur_worker], :] - z[id_scen, :]
        
        task = SubproblemTask(
            pb.scenarios[id_scen],
            id_scen,
            pb.build_subpb,
            μ,
            v,
        )
        
        ## send task to worker
        put!(work_channels[cur_worker], task)
        worker_to_launchit[cur_worker] = it
        worker_to_scen[cur_worker] = id_scen

        steplength = norm(step)
        # invariants, indicators, prints
        if mod(it, params[:print_step]) == 0
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            
            @printf "%3i    %.10e % .16e\n" it steplength objval
        end
        
        it += 1
    end
    
    ## Get final solution
    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    
    @printf "%3i    %.10e % .16e\n" it steplength objval
    println("Computation time: ", time() - t_init)
    
    ## Terminate all workers
    printstyled("\nTerminating nodes...\n", color=:red)
    terminate_workers(pb, work_channels, remotecalls_futures)

    return x_feas
end
