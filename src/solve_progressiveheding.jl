function subproblem_solve(pb, id_scen, u_scen, x_scen, μ, params)
    n = sum(length.(pb.stage_to_dim))
    
    ## Regalurized problem
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)
    
    # Augmented lagragian subproblem full objective
    obj += dot(u_scen, y) + (1/2*μ) * sum((y[i]-x_scen[i])^2 for i in 1:n)
    @objective(model, Min, obj)
    
    optimize!(model)

    y_new = JuMP.value.(y)
    return y_new
end

"""
nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})

Store in `x` the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
    @assert size(x) == size(y) "nonanticipatory_projection!(): input and output arrays should have same size"
    depth_to_partition = get_partitionbydepth(pb.scenariotree)

    for (stage, scen_partition) in enumerate(depth_to_partition)
        stage_dims = pb.stage_to_dim[stage]
        for scen_set in scen_partition
            averaged_traj = sum(pb.probas[i]*y[i, stage_dims] for i in scen_set) / sum(pb.probas[i] for i in scen_set)
            for scen in scen_set
                x[scen, stage_dims] = averaged_traj
            end
        end
    end
    return
end

"""
nonanticipatory_projection(pb::Problem, y::Matrix{Float64})

Compute the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection(pb::Problem, y::Matrix{Float64})
    x = zeros(size(y))
    nonanticipatory_projection!(x, pb, y)
    return x
end


"""
    solve_progressivehedging(pb::Problem)

Run the classical Progressive Hedging scheme on problem `pb`.
"""
function solve_progressivehedging(pb::Problem)
    println("--------------------------------------------------------")
    println("--- Progressive Hedging - sequential")
    println("--------------------------------------------------------")
    
    # parameters
    μ = 3
    params = Dict(
        :max_iter => 100,
        :print_step => 10,
        :ϵ_prim => 1e-3,
        :ϵ_dual => 1e-3,
    )

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    y = zeros(Float64, nscenarios, n)
    x = zeros(Float64, nscenarios, n)
    u = zeros(Float64, nscenarios, n)
    primres = dualres = Inf
    
    # Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    it = 0
    @printf " it   primal res        dual res            dot(x,u)   objective\n"
    while !(primres < params[:ϵ_prim] && dualres < params[:ϵ_dual]) && it < params[:max_iter]
        u_old = copy(u)

        # Subproblem solves
        for id_scen in 1:nscenarios
            y[id_scen, :] = subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, params)
        end

        # projection on non anticipatory subspace
        nonanticipatory_projection!(x, pb, y)

        # multiplier update
        u += (1/μ) * (y-x)

        
        # invariants, indicators
        primres = norm(pb, x-y)
        dualres = (1/μ) * norm(pb, u - u_old)
        if mod(it, params[:print_step]) == 0
            objval = objective_value(pb, x)
            dot_xu = dot(pb, x, u)
        
            @printf "%3i   %.10e  %.10e   % .3e % .16e\n" it primres dualres dot_xu objval
        end

        it += 1
    end

    ## Final print
    objval = objective_value(pb, x)
    dot_xu = dot(pb, x, u)

    @printf "%3i   %.10e  %.10e   % .3e % .16e\n" it primres dualres dot_xu objval

    return x
end