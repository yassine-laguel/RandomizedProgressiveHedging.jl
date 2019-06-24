using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.status()

using RPH, Ipopt
using Juniper, Cbc

include("build_hydrothermalscheduling.jl")
include("build_hydrothermalscheduling_extended.jl")
include("build_hydrothermalscheduling_milp.jl")

function main()
    nstages = 5
    # pb = build_hydrothermal_problem(nstages = nstages)
    pb = build_hydrothermalextendedmilp_problem(nstages = nstages, ndams=1)

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    # y_direct = solve_direct(pb)
    # println("\nDirect solve output is:")
    # display(y_direct)
    # println("")

    optimizer = Juniper.Optimizer
    optimizer_params = Dict{Symbol, Any}()
    optimizer_params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
    optimizer_params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)
    optimizer_params[:log_levels] = []

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    # y_PH = solve_progressivehedging(pb, maxtime=20, ϵ_primal=1e-4, printstep=1, optimizer=optimizer, optimizer_params=optimizer_params)
    # println("\nSequential solve output is:")
    # display(y_PH);
    
    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_sync = solve_randomized_sync(pb, maxtime=20, printstep=10, optimizer=optimizer, optimizer_params=optimizer_params)
    # println("\nSynchronous solve output is:")
    # display(y_sync);
    
    # #########################################################
    # ## Problem solve: asynchronous (parallelized) version of PH
    y_async = solve_randomized_async(pb, maxtime=20, printstep=100)
    # println("Asynchronous solve output is:")
    # display(y_async);

    @show norm(y_sync - y_PH)
    @show norm(y_sync - y_sync)
    @show norm(y_sync - y_async)

    return
end

main()
