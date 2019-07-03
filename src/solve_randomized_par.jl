"""
    SubproblemTask{T}

TODO
"""
struct SubproblemTask{T}
    taskid::Int
    scenario::T
    id_scenario::ScenarioId
    v_scen::Vector{Float64}
end

"""
    do_remote_work(inwork::RemoteChannel, outres::RemoteChannel)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`v_scen`)' where 'f_s' is the cost function associated with 
the scenario `id_scen`.
"""
function do_remote_work(inwork::RemoteChannel, outres::RemoteChannel, paramschan::RemoteChannel)
    params = take!(paramschan)
    μ = params[:μ]
    build_fs = params[:build_fs]

    while true
        try
            t0 = time()
            subpbtask::SubproblemTask = take!(inwork)

            if subpbtask.taskid == -1  # Work finished
                return
            end

            model = Model(with_optimizer(params[:optimizer]; params[:optimizer_params]...))
            
            # Get scenario objective function, build constraints in model
            y, obj, ctrref = build_fs(model, subpbtask.scenario, subpbtask.id_scenario)
            
            obj += (1/2*μ) * sum((y[i] - subpbtask.v_scen[i])^2 for i in 1:length(y))
            
            @objective(model, Min, obj)

            optimize!(model)
            if (primal_status(model) !== MOI.FEASIBLE_POINT) || (dual_status(model) !== MOI.FEASIBLE_POINT)
                @warn "Subproblem of scenario $id_scen " primal_status(model) dual_status(model) termination_status(model)
            end

            put!(outres, (JuMP.value.(y), subpbtask.id_scenario))
            # println(time()-t0)
        catch e
            println("Worker error:")
            println(e)
        end
    end
end


"""
    init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario

TODO
"""
function init_workers(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario
    if workers() == Vector([1])
        @error "No workers available. Returning"
        return
    end

    nworkers = length(workers())

    work_channel = RemoteChannel(()->Channel{SubproblemTask{T}}(nworkers))
    results_channel = RemoteChannel(()->Channel{Tuple{Vector{Float64}, Int}}(nworkers))
    params_channel = RemoteChannel(()->Channel{Dict{Symbol,Any}}(nworkers))
    
    remotecalls_futures = OrderedDict(worker_id => remotecall(do_remote_work, 
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
    terminate_workers(pb, work_channel, remotecalls_futures)

TODO
"""
function terminate_workers(pb, work_channel, remotecalls_futures)
    closing_task = SubproblemTask(
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
    randomizedasync_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)

TODO
"""
function randomizedasync_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)
    printlev>0 && print("Initialisation... ")
    xz_scen = zeros(get_scenariodim(pb))
    cur_scen = 1
    ## First, feed as many scenarios as there are workers
    for id_worker in 1:min(length(workers()), pb.nscenarios)
        put!(work_channel, SubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
        cur_scen += 1
    end

    it = 0
    while it < pb.nscenarios
        ## Taking current result
        y, id_scen = take!(results_channel)
        z[id_scen, :] = y

        ## Submitting new task if necessary
        if cur_scen <= pb.nscenarios
            put!(work_channel, SubproblemTask(cur_scen, pb.scenarios[cur_scen], cur_scen, xz_scen))
            cur_scen += 1
        end
        it += 1
    end

    nonanticipatory_projection!(z, pb, z)
    
    printlev>0 && println("done")
    return it
end


function init_hist!(hist, printstep)
    !isnothing(hist) && (hist[:functionalvalue] = Float64[])
    !isnothing(hist) && (hist[:time] = Float64[])
    !isnothing(hist) && haskey(hist, :approxsol) && (hist[:dist_opt] = Float64[])
    !isnothing(hist) && (hist[:logstep] = printstep)
    return
end

function get_defaultsubpbparams(pb, optimizer, optimizer_params, μ)
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
function solve_randomized_par(pb::Problem{T}; μ::Float64 = 3.0,
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
    printlev>0 && println("--- Randomized Progressive Hedging - parallel")
    printlev>0 && println("--------------------------------------------------------")
    
    # Variables
    nworkers = length(workers())
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = get_scenariodim(pb)
    
    x = zeros(Float64, nworkers, n)
    z = zeros(Float64, nscenarios, n)
    step = zeros(Float64, n)
    steplength = Inf
    
    ## Random scenario sampling
    rng = MersenneTwister(isnothing(seed) ? 1234 : seed)
    scen_sampling_distrib = Categorical(isnothing(qdistr) ? pb.probas : qdistr)
    qmin = minimum(scen_sampling_distrib.p)      
    
    
    init_hist!(hist, printstep*min(nworkers,pb.nscenarios))
    
    subpbparams = get_defaultsubpbparams(pb, optimizer, optimizer_params, μ)

    ## Workers Initialization
    printlev>0 && println("Available workers: ", nworkers)
    work_channel, results_channel, params_channel, remotecalls_futures = init_workers(pb, subpbparams)
  
    it = 0
    tinit = time()
    printlev>0 && @printf "   it   residual            objective                \n"

    it = randomizedasync_initialization!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)
    x = copy(z)
    
    printlev>0 && @printf "%5i   %.10e   % .16e    \n" it 0.0 objective_value(pb, z)  
    



    while it < maxiter && time()-tinit < maxtime
        
        z_prev = copy(z)
        
        ## Draw a scenario, build v, send task
        tab_id_scen = zeros(Int, nworkers)
        for i in 1:nworkers
            id_scen = rand(rng, scen_sampling_distrib)
            while id_scen in tab_id_scen[1:i-1]
                id_scen = rand(rng, scen_sampling_distrib)
            end
            tab_id_scen[i] = id_scen
        end

        ## For all available workers
        for id_scen in tab_id_scen
            
            v = 2*x[id_scen, :] - z[id_scen, :]

            put!(work_channel, SubproblemTask(id_scen, pb.scenarios[id_scen], id_scen, v))
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

            !isnothing(hist) && push!(hist[:functionalvalue], objval)
            !isnothing(hist) && push!(hist[:time], time() - tinit)
            !isnothing(hist) && haskey(hist, :approxsol) && size(hist[:approxsol])==size(x) && push!(hist[:dist_opt], norm(hist[:approxsol] - x))                
        end
        
        it += 1
    end

    ## Get final solution

    objval = objective_value(pb, x)
    
    printlev>0 && mod(it, printstep) != 1 && @printf "%5i   %.10e   % .16e\n" it steplength objval 
    printlev>0 && println("Computation time: ", time() - tinit)
    
    ## Terminate all workers
    terminate_workers(pb, work_channel, remotecalls_futures)

    return x
end
