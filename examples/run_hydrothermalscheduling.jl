using Distributed
using RPH, Ipopt

include("build_hydrothermalscheduling_extended.jl")

function main()
    nstages = 5
    # pb = build_hydrothermal_problem(nstages = nstages)
    pb = build_hydrothermalextended_problem(nstages = nstages, ndams=5)

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pb, optimizer=GLPK.Optimizer)
    println("\nDirect solve output is:")
    # display(y_direct)
    # println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = solve_progressivehedging(pb, maxtime=20, ϵ_primal=1e-4, printstep=1)
    # println("\nSequential solve output is:")
    # display(y_PH);

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_sync = solve_randomized_sync(pb, maxtime=3, printstep=10)
    # println("\nSynchronous solve output is:")
    # display(y_sync);

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_par = solve_randomized_par(pb, maxtime=3, printstep=10)
    # println("\nSynchronous solve output is:")
    # display(y_sync);

    # #########################################################
    # ## Problem solve: asynchronous (parallelized) version of PH
    y_async = solve_randomized_async(pb, maxtime=15, printstep=100)
    # println("Asynchronous solve output is:")
    # display(y_async);

    @show norm(y_sync - y_PH)
    @show norm(y_sync - y_par)
    @show norm(y_sync - y_async)

    return
end

main()
