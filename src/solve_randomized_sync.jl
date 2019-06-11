# include("RPH.jl")
# include("PH_sequential.jl")

"""
get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)

Compute the average trajectory defined by scenario `id_scen` over strategy `z`.
"""
function get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)
    nstages = pb.nstages
    n = sum(length.(pb.stage_to_dim))

    averaged_traj = zeros(n)
    
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

    return averaged_traj
end

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
solve_randomized_sync(pb)

Run the Randomized Progressive Hedging scheme on problem `pb`.
"""
function solve_randomized_sync(pb)
    println("--------------------------------------------------------")
    println("--- Randomized Progressive Hedging - synchronous")
    println("--------------------------------------------------------")
    
    # parameters
    μ = 3
    params = Dict(
        :print_step=>5,
        :max_iter=>30,
    )

    scen_sampling_distrib = Categorical(pb.probas)

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    x = zeros(Float64, n)
    y = zeros(Float64, n)
    z = zeros(Float64, nscenarios, n)
    
    # Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    it = 0.0
    oldit = 0
    @printf " it   global residual   objective\n"
    while it < params[:max_iter] * nscenarios
        id_scen = rand(scen_sampling_distrib)

        ## Projection
        x = get_averagedtraj(pb, z, id_scen)

        ## Subproblem solve
        y = PH_sync_subpbsolve(pb, id_scen, 2*x-z[id_scen, :], μ, params)

        ## Global variable update
        z[id_scen, :] += (y - x)


        # invariants, indicators, prints
        if it > oldit + params[:print_step]

            # @btime x_feas = nonanticipatory_projection($pb, $z)
            x_feas = nonanticipatory_projection(pb, z)
            objval = objective_value(pb, x_feas)
            primres = norm(x-y)
            
            @printf "%3i   %.10e % .16e\n" it primres objval
            oldit += params[:print_step]
        end
        
        it += 1/nscenarios
    end

    x_feas = nonanticipatory_projection(pb, z)
    objval = objective_value(pb, x_feas)
    primres = norm(x-y)
    
    @printf "%3i   %.10e % .16e\n" it primres objval

    return x_feas
end