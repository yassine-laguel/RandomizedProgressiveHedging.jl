using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.status()

using RPH

include("build_hydrothermalscheduling.jl")

function main()
    nstages = 10
    # pb = build_hydrothermal_problem(nstages = nstages)
    pb = build_hydrothermal_problem_vscenario(nstages = nstages)

    println("Full problem is:")
    println(pb)
    
    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    # y_direct = solve_direct(pb, solver=with_optimizer(GLPK.Optimizer))
    # println("\nDirect solve output is:")
    # display(y_direct)
    # println("")

    # #########################################################
    # ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    # y_PH = solve_progressivehedging(pb)
    # println("\nSequential solve output is:")
    # display(y_PH)
    # println("")

    # #########################################################
    # ## Problem solve: synchronous (un parallelized) version of PH
    # y_synch = solve_randomized_sync(pb)
    # println("\nSynchronous solve output is:")
    # display(y_synch)
    
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
