"""
PH_sync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)

Solve and return the solution of the subproblem 'prox_(f_s/`μ`) (`xz_scen`)' where 'f_s' is the cost function associated with 
the scenario `id_scen`.
"""
function PH_sync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)
    n = sum(length.(pb.stage_to_dim))

    ## Regalurized problem
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

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
"""
function solve_randomized_sync(pb::Problem, kwargs...)
    println("--------------------------------------------------------")
    println("--- Randomized Progressive Hedging - synchronous")
    println("--------------------------------------------------------")
    
    # parameters
    μ = 3
    params = Dict(
        :print_step=>10,
        :max_iter=>400,
        :step_tol=>1e-4,
    )

    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))
    
    qproba = ones(nscenarios) / nscenarios
    scen_sampling_distrib = Categorical(pb.probas)


    # Optim Variables
    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    steplength = Inf
    
    # Initialization
    # z = subproblems per scenario

    it = 0
    @printf " it   global residual   objective\n"
    while !(steplength < params[:step_tol]) && it < params[:max_iter]
        id_scen = rand(scen_sampling_distrib)

        ## Projection
        get_averagedtraj!(x, pb, z, id_scen) #TODO: rename with proj ?

        ## Subproblem solve
        y = PH_sync_subpbsolve(pb, id_scen, 2*x-z[id_scen, :], μ, params)

        ## Global variable update
        z[id_scen, :] += (y - x)


        # invariants, indicators, prints
        if mod(it, params[:print_step]) == 0
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            steplength = norm(x-y)
            
            @printf "%3i   %.10e % .16e\n" it steplength objval
        end
        
        it += 1
    end

    ## Final print
    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    steplength = norm(x-y)
    
    @printf "%3i   %.10e % .16e\n" it steplength objval

    return x_feas
end
