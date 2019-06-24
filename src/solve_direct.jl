function build_global_objective!(model, pb::Problem, scen_to_obj, ::Type{RiskNeutral})
    @objective(model, Min, sum(proba_s * scen_to_obj[id_s] for (id_s, proba_s) in enumerate(pb.probas)))
    return nothing
end

function build_global_objective!(model, pb::Problem, scen_to_obj, cvar::CVar)
    nscenarios = pb.nscenarios

    @variable(model, η)
    @variable(model, maxvals[1:nscenarios])

    @objective(model, Min, η + sum(proba_s / (1-cvar.p) * maxvals[id_s] for (id_s, proba_s) in enumerate(pb.probas)))

    ## Constraint enforcing maxvals[i] = max(0, f_i(x_i)-η)
    @constraint(model, maxvals .>= 0)
    @constraint(model, [id_s=1:nscenarios], maxvals[id_s] >= scen_to_obj[id_s] - η )

    return nothing
end

"""
    x = solve_direct(pb::Problem; optimizer = GLPK.Optimizer, printlev=1)

Build the progressive hedging problem by explicitly laying out non-anticipatory 
constraints, and solve globally.

## Keyword arguments:
- `optimizer`: optimizer used for solve. Default is `Ipopt.Optimizer`.
- `printlev`: if 0, mutes output from the function (not solver). Default value is 1.
"""
function solve_direct(pb::Problem; riskmeasure = RiskNeutral, 
                                   optimizer = GLPK.Optimizer, 
                                   optimizer_params = Dict{Symbol, Any}(),
                                   printlev = 1)
    printlev>0 && println("--------------------------------------------------------")
    printlev>0 && println("--- Direct solve")
    printlev>0 && println("--------------------------------------------------------")

    model = Model(with_optimizer(optimizer; optimizer_params...))

    printlev>0 && println("Building global model...")
    ## Collect subproblems
    scen_to_vars = SortedDict()
    scen_to_obj = SortedDict()
    scen_to_ctr = SortedDict()

    for (id_scen, scen) in enumerate(pb.scenarios)
        y, obj, ctrref = pb.build_subpb(model, scen, id_scen)
        scen_to_vars[id_scen] = y
        scen_to_obj[id_scen] = obj
        scen_to_ctr[id_scen] = ctrref
    end

    # Layout global objective
    build_global_objective!(model, pb, scen_to_obj, riskmeasure)

    ## Non anticipatory constraints
    printlev>0 && println("Laying out nonanticipatory constraints...")
    depth_to_part = get_partitionbydepth(pb.scenariotree)
    for (depth, partition) in enumerate(depth_to_part)
        for subset in partition
            id_min, id_max = subset.start, subset.stop
            for scen_i in subset, scen_j in scen_i+1:id_max
                stage_subspace = pb.stage_to_dim[depth]
                @constraint(model, scen_to_vars[scen_i][stage_subspace] .- scen_to_vars[scen_j][stage_subspace] .== 0)
            end
        end
    end

    ## Optimization and solution
    printlev>0 && print("Optimization... ")
    optimize!(model)
    printlev>0 && println("Done.")

    y_sol = zeros(pb.nscenarios, sum(length.(pb.stage_to_dim)))
    for (id_scen, vars) in scen_to_vars
        y_sol[id_scen, :] = JuMP.value.(vars)
    end

    printlev>0 && @show termination_status(model)
    printlev>0 && @show primal_status(model)
    printlev>0 && @show dual_status(model)

    return y_sol
end
