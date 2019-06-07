using Ipopt, Distributed
include("../src/RPH.jl")
include("../src/testcases.jl")
include("../src/PH_direct.jl")
include("../src/PH_sequential.jl")
include("../src/PH_synchronous.jl")
include("../src/PH_asynchronous.jl")

@everywhere struct SimpleExScenario <: AbstractScenario
    trajcenter::Vector{Float64}
    constraintbound::Int
end


@everywhere function build_fs_Cs!(model::JuMP.Model, s::SimpleExScenario, id_scen::ScenarioId)
    n = length(s.trajcenter)
    
    y = @variable(model, [1:n], base_name="y_s$id_scen")
    
    objexpr = sum((y[i] - s.trajcenter[i])^2 for i in 1:n)
    
    ctrref = @constraint(model, y .<= s.constraintbound)
    
    return y, objexpr, ctrref
end


function main()

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
        [0.5, 0.25, 0.25],                  # scenario probabilities
        [1:1, 2:2, 3:3],                    # stage id to trajectory coordinates, required for projection
        stageid_to_scenpart                 # stage to scenario partition
    )

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = PH_direct_solve(pb)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = PH_sequential_solve(pb)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_synch = PH_synchronous_solve(pb)
    println("\nSynchronous solve output is:")
    display(y_synch)
    
    #########################################################
    ## Problem solve: asynchronous (parallelized) version of PH
    y_asynch = PH_asynchronous_solve(pb)
    println("ASynchronous solve output is:")
    display(y_asynch)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_synch)
    @show norm(y_direct - y_asynch)
    
    return
end

main()