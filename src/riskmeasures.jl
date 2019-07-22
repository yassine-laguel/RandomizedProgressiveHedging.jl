export cvar_problem

"""
    pb_cvar = cvar_problem(pb::Problem, cvar::CVar)

Build the problem with risk measure corresponding to `cvar`.
"""
function cvar_problem(pb::Problem, cvarlevel::Real)
    ## Define extended problem: change stage_to_dim
    new_dim_to_subspace = [1:2 for i in 1:pb.nstages]
    
    new_dim_to_subspace[1] = 1:pb.stage_to_dim[1].stop+1
    for stageid in 2:pb.nstages
        new_dim_to_subspace[stageid] = pb.stage_to_dim[stageid].start+1:pb.stage_to_dim[stageid].stop+1
    end

    function build_cvarproblem!(model::JuMP.Model, s::T, id_scen::RPH.ScenarioId) where T<:AbstractScenario
        y, objexpr, ctrref = pb.build_subpb(model, s, id_scen)

        η = @variable(model)
        maxval = @variable(model)

        cvarobj = η + maxval/(1-cvarlevel)

        con1 = @constraint(model, maxval >= 0)
        con2 = @constraint(model, maxval >= objexpr - η)
        return vcat(η, y), cvarobj, union(ctrref, [con1, con2])
    end


    ## lay out new objective function
    cvarpb = Problem(
        pb.scenarios,
        build_cvarproblem!,
        pb.probas,
        pb.nscenarios, 
        pb.nstages,
        new_dim_to_subspace,
        pb.scenariotree
    )

    return cvarpb
end