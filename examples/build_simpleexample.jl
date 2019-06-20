using Distributed
using DataStructures, LinearAlgebra, GLPK
using RPH, JuMP

@everywhere struct SimpleExScenario <: RPH.AbstractScenario
    trajcenter::Vector{Float64}
    constraintbound::Int
end

@everywhere function build_fs_Cs!(model::JuMP.Model, s::SimpleExScenario, id_scen::RPH.ScenarioId)
    n = length(s.trajcenter)

    y = @variable(model, [1:n], base_name="y_s$id_scen")
    objexpr = sum((y[i] - s.trajcenter[i])^2 for i in 1:n)
    ctrref = @constraint(model, y .<= s.constraintbound)
    
    return y, objexpr, ctrref
end

function build_simpleexample()
        #########################################################
    ## Problem definition

    scenario1 = SimpleExScenario([1, 1, 1], 3)
    scenario2 = SimpleExScenario([2, 2, 2], 3)
    scenario3 = SimpleExScenario([3, 3, 3], 3)

    # stage to scenario partition
    stageid_to_scenpart = [
        OrderedSet([BitSet(1:3)]),                      # Stage 1
        OrderedSet([BitSet(1), BitSet(2:3)]),           # Stage 2
        OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),  # Stage 3
    ]

    pb = Problem(
        [scenario1, scenario2, scenario3],  # scenarios array
        build_fs_Cs!,
        [0.5, 0.25, 0.25],                  # scenario probabilities
        [1:1, 2:2, 3:3],                    # stage id to trajectory coordinates, required for projection
        stageid_to_scenpart                 # stage to scenario partition
    )

    return pb
end