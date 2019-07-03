using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
# @everywhere Pkg.status()

using RPH, Ipopt
using Juniper, Cbc

include("build_hydrothermalscheduling.jl")
include("build_hydrothermalscheduling_extended.jl")
include("build_hydrothermalscheduling_milp.jl")

function main()
    pb = build_hydrothermalextendedmilp_problem(nstages = 2, ndams=2)

    println("Full problem is:")
    println(pb)

    optimizer = Juniper.Optimizer
    optimizer_params = Dict{Symbol, Any}()
    optimizer_params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
    optimizer_params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)
    optimizer_params[:log_levels] = []
    # optimizer_params[:log_levels] = [:Info, :Table]

    # #########################################################
    # ## Problem solve: build and solve complete problem, exponential in constraints
    # y_direct = solve_direct(pb, optimizer=optimizer, optimizer_params=optimizer_params)
    # println("\nDirect solve output is:")
    # display(y_direct)
    # println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = solve_progressivehedging(pb, maxtime=20, ϵ_primal=1e-8, ϵ_dual=1e-8, printstep=1, optimizer=optimizer, optimizer_params=optimizer_params)
    println("\nProgressive hedging solve output is:")
    display(y_PH);

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_sync = solve_randomized_sync(pb, maxtime=3, printstep=10, optimizer=optimizer, optimizer_params=optimizer_params)
    println("\nSynchronous solve output is:")
    display(y_sync);


    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_par = solve_randomized_par(pb, maxtime=20, printstep=10, optimizer=optimizer, optimizer_params=optimizer_params)
    println("\nParallel sollve output is:")
    display(y_par);


    #########################################################
    ## Problem solve: asynchronous (parallelized) version of PH
    y_async = solve_randomized_async(pb, maxtime=20, printstep=10, optimizer=optimizer, optimizer_params=optimizer_params)
    println("Asynchronous solve output is:")
    display(y_async);

    @show norm(y_sync - y_PH)
    @show norm(y_sync - y_sync)
    @show norm(y_sync - y_async)


    println("\nProgressive hedging solve output is:")
    display(y_PH);
    println("\nSynchronous solve output is:")
    display(y_sync);
    println("\nParallel sollve output is:")
    display(y_par);
    println("Asynchronous solve output is:")
    display(y_async);

    return
end

main()
