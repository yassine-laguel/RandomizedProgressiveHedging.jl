# ENV["OAR_NODEFILE"] = joinpath(".", "logdir", "config")
using Distributed, OarClusterManager

@assert basename(pwd())=="RPH.jl" "This script should be run from the RPH.jl folder."

GLOBAL_LOG_DIR = joinpath("/", "bettik", "PROJECTS", "pr-cvar", "RPH_num_exps")
# GLOBAL_LOG_DIR = joinpath(".", "logdir")

## Add all available workers
# !(workers() == Vector([1])) && (rmprocs(workers()); println("removing workers"))
# addprocs(get_ncoresmaster()-1)
# length(get_remotehosts())>0 && addprocs_oar(get_remotehosts())

## Load relevant packages in all workers
push!(LOAD_PATH, pwd())
using RPH, JuMP

using GLPK, Ipopt, LinearAlgebra
using DataStructures, Dates, DelimitedFiles, JSON

include("../examples/build_simpleexample.jl")
include("../examples/build_hydrothermalscheduling_extended.jl")
include("utils.jl")


function main()
    println("Entering main()")
    println("- Date is: ", String(Dates.format(now(), "yyyy_mm_dd-HHhMM")))
    println("  Path is: ", pwd())
    println("- ENV[\"OAR_NODEFILE\"]: ", get(ENV, "OAR_NODEFILE", ""))
    println("  # available cores on master node : ", get_ncoresmaster())
    println("  # available remote cores         : ", length(get_remotehosts()))

    ## Set up logging directory
    logdir = set_logdir()
    println("- Logging directory: ", logdir)
    println()


    ## Build problems to be solved
    problems = []

    # push!(problems, OrderedDict(
    #     :pbname => "simpleproblem",
    #     :pb => build_simpleexample(),
    # ))
    # push!(problems, OrderedDict(
    #     :pbname => "hydrothermal_7stages_20dams",
    #     :pb => build_hydrothermalextended_problem(;nstages=7, ndams=20),
    # ))
    push!(problems, OrderedDict(
        :pbname => "hydrothermal_12stages_5dams",
        :pb => build_hydrothermalextended_problem(;nstages=12, ndams=5),
    ))

    ## Set number of seeds to be tried
    maxtime = 3*60*60
    maxiter = 1e9
    seeds = 1:5

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
        :algoname => "randomized_sync",
        :fnsolve_symbol => :solve_randomized_sync,
        :maxtime => maxtime,
        :maxiter => maxiter,
        :seeds => seeds,
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
    runallalgs(;runasync=false)

    ## Solve
    for problem_descr in problems
        @show problem_descr
        pb = problem_descr[:pb]
        pbname = problem_descr[:pbname]
        println("\n- [", String(Dates.format(now(), "HHhMM SS")), "] Dealing with new problem:")
        println(pb)

        
        algo_to_seedhist = OrderedDict{String, Any}()


        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        println("[", String(Dates.format(now(), "HHhMM SS")), "] Long ph solve...")
        xsol = nothing
        fopt = nothing
        try
            xsol = solve_progressivehedging(pb, ϵ_primal=1e-10, ϵ_dual=1e-10, maxtime=4*60*60, maxiter=1e6, printstep=1)
            fopt = objective_value(pb, xsol)
        catch e
            println("Error during ph solve.")
            println(e)
        end


        ## Then, run other algs
        for algo_descr in algorithms
            println("\n  - [", String(Dates.format(now(), "HHhMM")), "] Running algo:")
            println(algo_descr[:fnsolve_symbol])
            println("    seeds:         ", algo_descr[:seeds])
            
            algo_to_seedhist[algo_descr[:algoname]] = OrderedDict()
            algo_to_seedhist[algo_descr[:algoname]][:fnsolve_symbol] = algo_descr[:fnsolve_symbol]
            algo_to_seedhist[algo_descr[:algoname]][:maxtime] = algo_descr[:maxtime]
            algo_to_seedhist[algo_descr[:algoname]][:maxiter] = algo_descr[:maxiter]
            algo_to_seedhist[algo_descr[:algoname]][:seeds] = algo_descr[:seeds]

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
