using Distributed

@assert basename(pwd())=="RPH.jl" "This script should be run from the RPH.jl folder."

GLOBAL_LOG_DIR = joinpath(".", "logdir")


## Load relevant packages in all workers
@everywhere push!(LOAD_PATH, pwd())
@everywhere using RPH, JuMP

using GLPK, Ipopt, LinearAlgebra
using DataStructures, Dates, DelimitedFiles, JSON

include("../examples/build_simpleexample.jl")
include("../examples/build_hydrothermalscheduling_extended.jl")
include("../examples/build_hydrothermalscheduling_milp.jl")
include("utils.jl")


function main()
    println("Entering main()")
    println("- Date is: ", String(Dates.format(now(), "yyyy_mm_dd-HHhMM")))
    println("  Path is: ", pwd())

    ## Set up logging directory
    logdir = set_logdir()
    println("- Logging directory: ", logdir)
    println()


    ## Build problems to be solved
    problems = []

    push!(problems, OrderedDict(
        :pbname => "simpleproblem",
        :pb => build_simpleexample(),
    ))
    push!(problems, OrderedDict(
        :pbname => "hydrothermal_5stages_20dams",
        :pb => build_hydrothermalextended_problem(;nstages=6, ndams=10),
    ))
    push!(problems, OrderedDict(
        :pbname => "hydrothermal_5stages_20dams",
        :pb => build_hydrothermalextendedmilp_problem(;nstages=6, ndams=10),
    ))


    nworkers = length(workers())
    nscenarios = 2^5

    ## Set number of seeds to be tried
    maxtime = 90
    maxiter = 400
    seeds = 1:4

    ## Build algorithms & params used for solve
    algorithms = []

    push!(algorithms, OrderedDict(
        :algoname => "progressivehedging",
        :fnsolve_symbol => :solve_progressivehedging,
        :maxtime => maxtime,
        :maxiter => maxiter,
        :seeds => 1,
        :printstep => 1,
    ))
    push!(algorithms, OrderedDict(
        :algoname => "randomized_sync0.001",
        :fnsolve_symbol => :solve_randomized_sync,
        :maxtime => maxtime,
        :maxiter => maxiter*nscenarios,
        :seeds => seeds,
        :mu => 0.001,
        :printstep => 20,
    ))
    push!(algorithms, OrderedDict(
        :algoname => "randomized_sync1",
        :fnsolve_symbol => :solve_randomized_sync,
        :maxtime => maxtime,
        :maxiter => maxiter*nscenarios,
        :seeds => seeds,
        :mu => 1,
        :printstep => 20,
    ))


    push!(algorithms, OrderedDict(
        :algoname => "randomized_par",
        :fnsolve_symbol => :solve_randomized_par,
        :maxtime => maxtime,
        :maxiter => maxiter*nscenarios/nworkers,
        :seeds => seeds,
        :printstep => 20,
    ))
    push!(algorithms, OrderedDict(
        :algoname => "randomized_async1",
        :fnsolve_symbol => :solve_randomized_async,
        :maxtime => maxtime,
        :maxiter => maxiter*nscenarios,
        :seeds => seeds,
        :mu => 1,
        :printstep => 20,
    ))

    push!(algorithms, OrderedDict(
        :algoname => "randomized_async0.0001",
        :fnsolve_symbol => :solve_randomized_async,
        :maxtime => maxtime,
        :maxiter => maxiter*nscenarios,
        :seeds => seeds,
        :mu => 0.0001,
        :printstep => 20,
    ))



    println("Experiment summarys:")
    println("  #problems:  ", length(problems))
    println("  algorithms: ", [a[:fnsolve_symbol] for a in algorithms])
    println("  seeds:      ", seeds)
    println()

    ## Logging object
    problem_to_algo = OrderedDict{String, Any}()
    problem_to_algo["problem_names"] = [pb[:pbname] for pb in problems]

    ## Run all algorithms once to precompile everything
    println("[", String(Dates.format(now(), "HHhMM SS")), "] Running algs once to precompile...")
    runallalgs()

    ## Solve
    for problem_descr in problems
        @show problem_descr
        pb = problem_descr[:pb]
        pbname = problem_descr[:pbname]
        println("\n- [", String(Dates.format(now(), "HHhMM SS")), "] Dealing with new problem:")
        println(pb)

        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        println("[", String(Dates.format(now(), "HHhMM SS")), "] Long ph solve...")
        xsol = nothing
        fopt = nothing
        try
            xsol = solve_progressivehedging(pb, ϵ_primal=1e-10, ϵ_dual=1e-10, maxtime=2*maxtime, maxiter=2*maxiter, printstep=10)
            fopt = objective_value(pb, xsol)
        catch e
            println("Error during ph solve.")
            println(e)
        end


        ## Then, run all algs
        algo_to_seedhist = OrderedDict{String, Any}()

        for algo_descr in algorithms
            println("\n  - [", String(Dates.format(now(), "HHhMM")), "] Running algo:")
            println(algo_descr[:fnsolve_symbol])
            println("    seeds:         ", algo_descr[:seeds])
            
            algo_to_seedhist[algo_descr[:algoname]] = OrderedDict()
            algo_to_seedhist[algo_descr[:algoname]][:fnsolve_symbol] = algo_descr[:fnsolve_symbol]
            algo_to_seedhist[algo_descr[:algoname]][:maxtime] = algo_descr[:maxtime]
            algo_to_seedhist[algo_descr[:algoname]][:maxiter] = algo_descr[:maxiter]
            algo_to_seedhist[algo_descr[:algoname]][:seeds] = algo_descr[:seeds]
            algo_to_seedhist[algo_descr[:algoname]][:nworkers] = length(workers())

            fnsolve = eval(algo_descr[:fnsolve_symbol])

            for seed in algo_descr[:seeds]
                println("  - [", String(Dates.format(now(), "HHhMM SS")), "] Solving for seed $seed")

                ## Set up log object
                hist = OrderedDict{Symbol, Any}()
                !isnothing(fopt) && (hist[:fopt] = fopt)
                !isnothing(xsol) && (hist[:approxsol] = xsol)

                println("\n--------------------------------------------------------")
                fnsolve(pb; maxtime = algo_descr[:maxtime], 
                            maxiter = algo_descr[:maxiter], 
                            printstep = algo_descr[:printstep],
                            hist = hist, 
                            seed = seed, 
                            ϵ_primal = 1e-10, 
                            ϵ_dual = 1e-10)
                println("--------------------------------------------------------\n")

                # Log seed, algo
                haskey(hist, :approxsol) && delete!(hist, :approxsol)
                algo_to_seedhist[algo_descr[:algoname]][seed] = hist
            end
        end

        problem_to_algo[pbname] = algo_to_seedhist
        problem_to_algo["nstages"] = pb.nstages
        problem_to_algo["nscenarios"] = pb.nscenarios
        problem_to_algo["subpbdim"] = sum(length.(pb.stage_to_dim))
        problem_to_algo["nworkers"] = length(workers())
    end

    ## Write log information
    println("[", String(Dates.format(now(), "HHhMM SS")), "] All computations are completed. Writing logs.")
    writelogs(problem_to_algo, logdir)

    println("[", String(Dates.format(now(), "HHhMM SS")), "] All done.")
    return
end

main()
