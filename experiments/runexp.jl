# ENV["OAR_NODEFILE"] = joinpath(".", "logdir", "config")
using Distributed, OarClusterManager

@assert basename(pwd())=="RPH.jl" "This script should be run from the RPH.jl folder."
@assert get_ncoresmaster()+length(get_remotehosts())>1 "At least one worker should ba available for async alg."

GLOBAL_LOG_DIR = joinpath("/", "bettik", "PROJECTS", "pr-cvar", "RPH_num_exps")
# GLOBAL_LOG_DIR = joinpath(".", "logdir")

## Add all available workers
!(workers() == Vector([1])) && (rmprocs(workers()); println("removing workers"))
addprocs(get_ncoresmaster()-1)
length(get_remotehosts())>0 && addprocs_oar(get_remotehosts())

## Load relevant packages in all workers
@everywhere push!(LOAD_PATH, pwd())
@everywhere using RPH, JuMP

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
    push!(problems, OrderedDict(
        :pbname => "hydrothermal_10stages_20dams",
        :pb => build_hydrothermalextended_problem(;nstages=10, ndams=20),
    ))

    ## Build algorithms & params used for solve
    algorithms = []

    push!(algorithms, OrderedDict(
        :algoname => "progressivehedging",
        :fnsolve_symbol => :solve_progressivehedging,
        :maxtime => 60*60,
        :maxiter => 1e9,
    ))
    push!(algorithms, OrderedDict(
        :algoname => "randomized_sync",
        :fnsolve_symbol => :solve_randomized_sync,
        :maxtime => 60*60,
        :maxiter => 1e9,
    ))
    push!(algorithms, OrderedDict(
        :algoname => "randomized_async",
        :fnsolve_symbol => :solve_randomized_async,
        :maxtime => 60*60,
        :maxiter => 1e9,
    ))

    ## Set number of seeds to be tried
    seeds = 1:5

    println("Experiment summarys:")
    println("  #problems:  ", length(problems))
    println("  algorithms: ", [a[:fnsolve_symbol] for a in algorithms])
    println("  seeds:      ", seeds)
    println()

    ## Logging object
    problem_to_algo = OrderedDict{String, Any}()
    problem_to_algo["problem_names"] = [pb[:pbname] for pb in problems]

    ## Run all algorithms once to precompile everything
    println("Running algs once to precompile...")
    runallalgs()

    ## Solve
    for problem_descr in problems
        @show problem_descr
        pb = problem_descr[:pb]
        pbname = problem_descr[:pbname]
        println("\n- [", String(Dates.format(now(), "HHhMM")), "] Dealing with new problem:")
        println(pb)

        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        println("Long solve for solution...")
        xsol = nothing
        fopt = nothing
        try
            xsol = solve_progressivehedging(pb, ϵ_primal=1e-5, ϵ_dual=1e-5, maxtime=3*60*60, maxiter=1e6, printstep=100)
            fopt = objective_value(pb, xsol)
        catch e
            println("Error during ph solve.")
            println(e)
        end


        ## Then, run all algs
        algo_to_seedhist = OrderedDict{String, Any}()

        for algo_descr in algorithms
            println("  - [", String(Dates.format(now(), "HHhMM")), "] Running algo:")
            println(algo_descr[:fnsolve_symbol])
            println("    nseeds:         ", seeds)
            
            algo_to_seedhist[algo_descr[:algoname]] = OrderedDict()
            algo_to_seedhist["seeds"] = collect(seeds)
            
            for seed in seeds
                println("  - Solving for seed $seed")

                ## Set up log object
                hist = OrderedDict{Symbol, Any}()
                !isnothing(fopt) && (hist[:fopt] = fopt)
                !isnothing(fopt) && (hist[:approxsol] = xsol)

                fnsolve = eval(algo_descr[:fnsolve_symbol])
                println("\n--------------------------------------------------------")
                fnsolve(pb; maxtime=algo_descr[:maxtime], maxiter=algo_descr[:maxiter], hist=hist, seed=seed)
                println("--------------------------------------------------------\n")

                # Log seed, algo
                haskey(hist, :approxsol) && delete!(hist, :approxsol)
                algo_to_seedhist[algo_descr[:algoname]][seed] = hist
                algo_to_seedhist[algo_descr[:algoname]][:fnsolve_symbol] = algo_descr[:fnsolve_symbol]
                algo_to_seedhist[algo_descr[:algoname]][:maxtime] = algo_descr[:maxtime]
                algo_to_seedhist[algo_descr[:algoname]][:maxiter] = algo_descr[:maxiter]
            end
        end

        problem_to_algo[pbname] = algo_to_seedhist
        problem_to_algo["nstages"] = pb.nstages
        problem_to_algo["nscenarios"] = pb.nscenarios
        problem_to_algo["subpbdim"] = sum(length.(pb.stage_to_dim))
        problem_to_algo["nworkers"] = length(workers())
    end

    ## Write log information
    println("All computations are completed. Writing logs.")
    writelogs(problem_to_algo, logdir)

    println("All done.")
    return
end

main()
