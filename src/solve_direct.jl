"""
    x = solve_direct(pb::Problem; optimizer = GLPK.Optimizer, printlev=1)

Build the progressive hedging problem by explicitly laying out non-anticipatory
constraints, and solve globally.

## Keyword arguments:
- `optimizer`: optimizer used for solve. Default is `GLPK.Optimizer`.
- `optimizer_params`: a `Dict{Symbol, Any}` storing parameters for the optimizer.
- `printlev`: if 0, mutes output from the function (not solver). Default value is 1.
"""
function solve_direct(pb::Problem; optimizer = GLPK.Optimizer,
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
    @objective(model, Min, sum(proba_s * scen_to_obj[id_s] for (id_s, proba_s) in enumerate(pb.probas)))

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
