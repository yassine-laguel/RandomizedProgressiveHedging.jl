using Distributed

@everywhere using JuMP, RPH

include("build_simpleexample.jl")
using Ipopt

function main()
    pb = build_simpleexample()

    hist=OrderedDict{Symbol, Any}(
        :approxsol => [1.75  1.0  1.0
        1.75  2.5  2.0
        1.75  2.5  2.99995]
    )

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pb, optimizer = Ipopt.Optimizer, printlev=0)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = solve_progressivehedging(pb, maxtime=1, printstep=1, hist=hist)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_sync = solve_randomized_sync(pb, maxtime=1.5, printstep=10, hist=hist)
    println("\nSynchronous solve output is:")
    display(y_sync)

    #########################################################
    ## Problem solve: asynchronous (parallelized) version of PH
    y_par = solve_randomized_par(pb, maxtime=1.5, printstep=10, hist=hist)
    println("\nRandom Par solve output is:")
    display(y_par)


    #########################################################
    ## Problem solve: asynchronous (parallelized) version of PH
    y_async = solve_randomized_async(pb, maxtime=1.5, printstep=10, hist=hist)
    println("Asynchronous solve output is:")
    display(y_async)

    return
end

main()
