include("RPH.jl")

using LinearAlgebra

@everywhere struct MyScenario <: AbstractScenario
    trajcenter::Vector{Float64}
end

@everywhere struct MySecondScenario <: AbstractScenario
    trajcenter::Vector{Float64}
    constraintbound::Float64
end


@everywhere function build_fs_Cs!(model::JuMP.Model, s::MyScenario, id_scen::ScenarioId)
    n = length(s.trajcenter)
    y = @variable(model, [1:n], base_name="y_s$id_scen")
    
    objexpr = sum((y[i] - s.trajcenter[i])^2 for i in 1:n)

    return y, objexpr, nothing
end

@everywhere function build_fs_Cs!(model::JuMP.Model, s::MySecondScenario, id_scen::ScenarioId)
    n = length(s.trajcenter)
    y = @variable(model, [1:3], base_name="y_s$id_scen")
    
    objexpr = sum((y[i] - s.trajcenter[i])^2 for i in 1:n)
    
    ctrref = @constraint(model, y .<= s.constraintbound)
    
    return y, objexpr, ctrref
end


function makeproblem()
    n = 3

    ## Scenario 1
    scenario1 = MyScenario([1, 1, 1])

    ## Scenario 2
    scenario2 = MyScenario([2, 2, 2])

    ## Scenario 3
    scenario3 = MySecondScenario([3, 3, 3], 4)

    stageid_to_scenpart = [
        OrderedSet([BitSet(1:3)]),
        OrderedSet([BitSet(1), BitSet(2:3)]),
        OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),
    ]


    pb = Problem(
        [scenario1, scenario2, scenario3],
        [0.5, 0.25, 0.25],
        [1:1, 2:2, 3:3],
        stageid_to_scenpart
    )

    return pb
end