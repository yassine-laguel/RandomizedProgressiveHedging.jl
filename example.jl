using Ipopt, Distributed
include("src/RPH.jl")
include("src/testcases.jl")
include("src/PH_direct.jl")
include("src/PH_sequential.jl")
include("src/PH_synchronous.jl")

@everywhere struct MyScenario <: AbstractScenario
    trajcenter::Vector{Float64}
end


@everywhere function build_fs_Cs!(model::JuMP.Model, s::MyScenario, id_scen::ScenarioId)
    n = length(s.trajcenter)
    
    y = @variable(model, [1:n], base_name="y_s$id_scen")
    
    objexpr = sum((y[i] - s.trajcenter[i])^2 for i in 1:n)
    
    ctrref = @constraint(model, y .<= s.constraintbound)
    
    return y, objexpr, ctrref
end


function main()

    #########################################################
    ## Problem definition

    scenario1 = MyScenario([1, 1, 1])
    scenario2 = MyScenario([2, 2, 2])
    scenario3 = MyScenario([3, 3, 3])

    # stage to scenario partition
    stageid_to_scenpart = [
        OrderedSet([BitSet(1:3)]),                      # Stage 1
        OrderedSet([BitSet(1), BitSet(2:3)]),           # Stage 2
        OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),  # Stage 3
    ]

    pb = Problem(
        [scenario1, scenario2, scenario3],  # scenarios array
        [0.5, 0.25, 0.25],                  # scenario probabilities
        [1:1, 2:2, 3:3],                    # stage id to trajectory coordinates, required for projection
        stageid_to_scenpart                 # stage to scenario partition
    )

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_sol = PH_direct_solve(pb)
    println("\nDirect solve output is:")
    display(y_sol)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_sol = PH_sequential_solve(pb)
    println("\nSequential solve output is:")
    display(y_sol)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_sol = PH_synchronous_solve(pb)
    println("\nSynchronous solve output is:")
    display(y_sol)

    return
end

main()