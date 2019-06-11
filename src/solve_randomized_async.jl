# include("RPH.jl")

"""
get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)

Compute the average trajectory defined by scenario `id_scen` over strategy `z`.
"""
function get_averagedtraj(pb::Problem, z::Matrix, id_scen::ScenarioId)
    nstages = pb.nstages
    n = sum(length.(pb.stage_to_dim))
    
    averaged_traj = zeros(n)
    
    try
    scentree = pb.scenariotree
    stage = 1
    id_curnode = scentree.idrootnode

    while stage <= scentree.depth
        ## Get scenarios, dimension for current stage
        scen_set = scentree.vecnodes[id_curnode].scenarioset
        stage_dims = pb.stage_to_dim[stage]

        ## update x section with average
        averaged_traj[stage_dims] = sum(pb.probas[i] * z[i, stage_dims] for i in scen_set)
        averaged_traj[stage_dims] /= sum(pb.probas[i] for i in scen_set)

        ## find node of following stage
        stage += 1
        id_nextnode = nothing
        for id_child in scentree.vecnodes[id_curnode].childs
            if id_scen in scentree.vecnodes[id_child].scenarioset
                id_nextnode = id_child
                break
            end
        end
        isnothing(id_nextnode) && stage < scentree.depth && @error "Tree dive, node of depth $stage not found."
        id_curnode = id_nextnode
    end

    catch e
        println(e)
    end

    return averaged_traj
end




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
            # println(round(time()-t0, sigdigits=2), "\tWaiting for work...")
            subpbtask::SubproblemTask = take!(inwork)
            # println(round(time()-t0, sigdigits=2), "\t... Got job rel. to scenario $(subpbtask.id_scenario).")

            if subpbtask.id_scenario == -1
                # Work finished
                # println(round(time()-t0, sigdigits=3), "\tI am done!")
                return
            end

            # do work
            model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
            
            # Get scenario objective function, build constraints in model
            y, obj, ctrref = subpbtask.build_subpb(model, subpbtask.scenario, subpbtask.id_scenario)
            
            obj += (1/2*subpbtask.μ) * sum((y[i] - subpbtask.v_scen[i])^2 for i in 1:length(y))
            
            @objective(model, Min, obj)

            sleep(0.05)
        
            optimize!(model)

            put!(outres, JuMP.value.(y))        
            # println(round(time()-t0, sigdigits=3), "\tI finished work")
        catch e
            println("Worker error:")
            println(e)
        end
    end
end



"""
solve_randomized_sync(pb)

Run the Randomized Progressive Hedging scheme on problem `pb`.
"""
function solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario
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
    )
    
    ## need workers() -> 1:nworkers map for x matrix
    worker_to_wid = SortedDict{Int, Int}(worker=>wid for (wid, worker) in enumerate(workers()))

    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))
    
    qproba = ones(nscenarios) / nscenarios
    scen_sampling_distrib = Categorical(qproba)

    c = 0.9
    qmin = minimum(qproba)      
    τ = ceil(Int, nworkers*1.1)                       # Upper bound on delay
    η = c*nscenarios*qmin / (2*τ*sqrt(qmin) + 1)
    
    # Variables
    x = zeros(Float64, nworkers, n)
    y = zeros(Float64, n)
    v = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)

    t_init = time()
    

    # Variables Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    ## Workers Initialization
    work_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{SubproblemTask{T}}(3), worker_id) for worker_id in workers())
    results_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Vector{Float64}}(3), worker_id) for worker_id in workers())
    
    worker_to_scen = OrderedDict{Int, ScenarioId}()
    worker_to_launchit = OrderedDict{Int, ScenarioId}()
    
    printstyled("Launching remotecalls...\n", color=:red)
    remotecalls_futures = OrderedDict(worker_id => remotecall(do_remote_work, worker_id, work_channels[worker_id], results_channels[worker_id]) for worker_id in workers())
    printstyled("Done.\n", color=:red)

    ## Feeding every worker with one task
    for w_id in workers()
        ## Draw random scenario, build subproblem task
        id_scen = rand(scen_sampling_distrib)

        task = SubproblemTask(
            pb.scenarios[id_scen],
            id_scen,
            pb.build_subpb,
            μ,
            v,
        )

        worker_to_scen[w_id] = id_scen
        worker_to_launchit[w_id] = 0

        put!(work_channels[w_id], task)
    end


    it = 0
    @printf "\n it    global residual   objective\n"
    while it < params[:max_iter]

        ## Wait for a worker to complete its job
        ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
        while length(ready_workers) == 0
            ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
        end
        
        ## Get oldest worker id TODO
        # worker_to_launchit[w_id] = 0
        cur_worker = first(ready_workers)
        
        ## Collect result
        y = take!(results_channels[cur_worker])
        id_scen = worker_to_scen[cur_worker]

        ## z update
        step = 2 * η / (nscenarios * pb.probas[id_scen]) * (y - x[worker_to_wid[cur_worker], :])
        z[id_scen, :] += step
        
        ## Draw new scenario for worker, build v and task
        id_scen = rand(scen_sampling_distrib)
        
        x[worker_to_wid[cur_worker], :] = get_averagedtraj(pb, z, id_scen)
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

        
        # invariants, indicators, prints
        if mod(it, params[:print_step]) == 0
            
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            primres = norm(step)
            
            @printf "%3i    %.10e % .16e\n" it primres objval
        end
        
        it += 1
    end

    ## Terminate all workers
    printstyled("\nTerminating nodes...\n", color=:red)

    closing_task = SubproblemTask(
        pb.scenarios[1],
        -1,                 ## Stop signal
        pb.build_subpb,
        μ,
        v,
    )
    for (w_id, worker_wkchan) in work_channels
        put!(worker_wkchan, closing_task)
    end

    for (w_id, worker) in remotecalls_futures
        wait(worker)
    end


    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    primres = 1e99

    println("Computation time: ", time() - t_init)
    
    @printf "%3i   %.10e % .16e\n" it primres objval

    return x_feas
end
