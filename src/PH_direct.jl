include("RPH.jl")
using JuMP, GLPK

"""
PH_direct_solve(pb::Problem)

Build the progressive hedging problem by explicitly laying out
non-anticipatory constraints, and solve globally. 
"""
function PH_direct_solve(pb::Problem)
    model = Model(with_optimizer(GLPK.Optimizer))

    ## Collect subproblems
    scen_to_vars = SortedDict()
    scen_to_obj = SortedDict()
    scen_to_ctr = SortedDict()

    for (id_scen, scen) in enumerate(pb.scenarios)
        y, obj, ctrref = build_fs_Cs!(model, scen, id_scen)
        scen_to_vars[id_scen] = y
        scen_to_obj[id_scen] = obj
        scen_to_ctr[id_scen] = ctrref
    end

    # Layout global objective
    @objective(model, Min, sum(proba_s * scen_to_obj[id_s] for (id_s, proba_s) in enumerate(pb.probas)))

    ## Non anticipatory constraints
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
    optimize!(model)
    
    y_sol = zeros(pb.nscenarios, sum(length.(pb.stage_to_dim)))
    for (id_scen, vars) in scen_to_vars
        y_sol[id_scen, :] = JuMP.value.(vars)
    end

    return y_sol
end
