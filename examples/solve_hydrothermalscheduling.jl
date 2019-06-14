include("build_hydrothermalscheduling.jl")

function main()
    nstages = 5
    # pb = build_hydrothermal_problem(nstages = nstages)
    pb = build_hydrothermal_problem_vscenario(nstages = nstages)

    println("Full problem is:")
    println(pb)
    # @show Base.summarysize(pb)
    # @show Base.summarysize(pb.scenarios)
    # @show Base.summarysize(pb.build_subpb)
    # @show Base.summarysize(pb.probas)
    # @show Base.summarysize(pb.nscenarios)
    # @show Base.summarysize(pb.nstages)
    # @show Base.summarysize(pb.stage_to_dim)
    # @show Base.summarysize(pb.scenariotree)


    # println()
    # stree = pb.scenariotree

    # @show Base.summarysize(stree)
    # @show Base.summarysize(stree.idrootnode)
    # @show Base.summarysize(stree.depth)
    # @show Base.summarysize(stree.vecnodes)
    # @show length(stree.vecnodes)
    # println()
    # @show Base.summarysize(stree.vecnodes[1])
    # @show Base.summarysize(stree.vecnodes[2])
    # @show Base.summarysize(stree.vecnodes[end])

    # println()
    # node = stree.vecnodes[end]
    # @show Base.summarysize(node)
    # @show Base.summarysize(node.father)
    # @show Base.summarysize(node.childs)
    # @show Base.summarysize(node.scenarioset)

    # node = stree.vecnodes[3]
    # @show Base.summarysize(node)
    # @show Base.summarysize(node.father)
    # @show Base.summarysize(node.childs)
    # @show node.childs
    # @show Base.summarysize(convert(Array{Int64}, node.childs))
    # @show Base.summarysize(node.scenarioset)
    
    
    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pb, solver=with_optimizer(GLPK.Optimizer))
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = solve_progressivehedging(pb)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_synch = solve_randomized_sync(pb)
    println("\nSynchronous solve output is:")
    display(y_synch)
    
    #########################################################
    ## Problem solve: asynchronous (parallelized) version of PH

    y_asynch = solve_randomized_async(pb)
    println("ASynchronous solve output is:")
    display(y_asynch)

    # @show norm(y_direct - y_PH)
    # @show norm(y_direct - y_synch)
    # @show norm(y_direct - y_asynch)

    return
end

main()
