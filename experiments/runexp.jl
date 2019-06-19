using Distributed, OarClusterManager

# GLOBAL_LOG_DIR = joinpath("/", "bettik", "PROJECTS", "pr-cvar")
GLOBAL_LOG_DIR = joinpath(".", "logdir")
ENV["OAR_NODEFILE"] = joinpath(".", "logdir", "config")

## Add all available workers
!(workers() == Vector([1])) && (rmprocs(workers()); println("removing workers"))
addprocs(get_ncoresmaster()-1)
length(get_remotehosts())>0 && addprocs_oar(get_remotehosts())

## Load relevant packages in all workers
@everywhere push!(LOAD_PATH, pwd())
@everywhere println(LOAD_PATH)
@everywhere using RPH, JuMP

using DataStructures, GLPK, LinearAlgebra
using Dates, DelimitedFiles, JSON

include("build_simpleexemple.jl")
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

    push!(problems, (
        pbname = "simpleproblem",
        pb = build_simpleexample(),
    ))

    ## Build algorithms & params used for solve
    algorithms = []

    push!(algorithms, (
        algoname = "solvedirect",
        fnsolve_symbol = :solve_direct,
        tlim = nothing,
        itmax = nothing,
    ))
    push!(algorithms, (
        algoname = "progressivehedging",
        fnsolve_symbol = :solve_progressivehedging,
        tlim = nothing,
        itmax = nothing,
    ))
    push!(algorithms, (
        algoname = "progressivehedging",
        fnsolve_symbol = :solve_randomized_sync,
        tlim = nothing,
        itmax = nothing,
    ))
    push!(algorithms, (
        algoname = "progressivehedging",
        fnsolve_symbol = :solve_randomized_async,
        tlim = nothing,
        itmax = nothing,
    ))

    ## Set number of seeds to be tried
    seeds = 1:2

    println("#problems:  ", length(algorithms))
    println("algorithms: ", [a.fnsolve_symbol for a in algorithms])
    println("seeds:      ", seeds)
    println()

    ## Logging object
    problem_to_algo = OrderedDict{Any, OrderedDict}()

    ## Run all algorithms once to precompile everything
    runallalgs()

    ## Solve
    for problem_descr in problems
        pb = problem_descr.pb
        println("- [", String(Dates.format(now(), "HHhMM")), "] Dealing with new problem:")
        println(pb)

        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        xsol, fopt = 0, 0

        algo_to_seedhist = OrderedDict{Any, OrderedDict}()

        ## Run algorithms x seed
        for algo_descr in algorithms
            println("- [", String(Dates.format(now(), "HHhMM")), "] Running algo:")
            println(algo_descr.fnsolve_symbol)
            println("nseeds:         ", seeds)
            
            algo_to_seedhist[algo_descr] = OrderedDict()
            
            for seed in seeds
                println("  - Solving for seed $seed")
                hist = OrderedDict{Symbol, Any}()
                
                fnsolve = eval(algo_descr.fnsolve_symbol)
                println("\n--------------------------------------------")
                fnsolve(pb; tlim=algo_descr.tlim, itmax=algo_descr.itmax, hist=hist, xsol=xsol, fopt=fopt, seed=seed)
                println("--------------------------------------------\n")

                # Log seed, algo
                algo_to_seedhist[algo_descr][seed] = hist
            end
        end

        problem_to_algo[problem_descr] = algo_to_seedhist
    end

    ## Write log information
    print("All computations are completed. Writing logs.")
    writelogs(problem_to_algo, logdir)

    println("All done.")
    return
end

main()
