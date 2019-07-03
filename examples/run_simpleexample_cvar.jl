using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
# @everywhere Pkg.status()

@everywhere using JuMP, RPH

include("build_simpleexample.jl")
using RPH, Ipopt

function main()
    pb = build_simpleexample()

    hist=OrderedDict{Symbol, Any}(
        :approxsol =>   [0.468975  1.72212  1.0      1.0    
                         0.468975  1.72212  2.62547  2.0    
                         0.468975  1.72212  2.62547  3.0]
    )

    println("Full problem is:")
    println(pb)

    cvar = CVar(0.1)
    pbcvar = cvar_problem(pb, cvar)
    
    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pbcvar, optimizer = Ipopt.Optimizer, printlev=0)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = solve_progressivehedging(pbcvar, maxiter=150, maxtime=40, printstep=5, hist=hist)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_synch = solve_randomized_sync(pbcvar, maxiter=600, maxtime=40, printstep=100)
    println("\nSynchronous solve output is:")
    display(y_synch)
    
    # #########################################################
    # ## Problem solve: asynchronous (parallelized) version of PH
    y_asynch = solve_randomized_async(pbcvar, maxiter=1e5, maxtime=20, printstep=100, hist=hist)
    println("Asynchronous solve output is:")
    display(y_asynch)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_synch)
    @show norm(y_direct - y_asynch)
    
    return
end

main()