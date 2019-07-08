using JSON, Dates, DataStructures
using Distributed, OarClusterManager
using JuMP

include("utils.jl")

function execute_algs_on_problems(problems, algorithms)
    println("Entering execute_algs_on_problems()")
    println("- Date is: ", String(Dates.format(now(), "yyyy_mm_dd-HHhMM")))
    println("  Path is: ", pwd())
    println("- ENV[\"OAR_NODEFILE\"]: ", get(ENV, "OAR_NODEFILE", ""))
    println("  # available cores on master node : ", get_ncoresmaster())
    println("  # available remote cores         : ", length(get_remotehosts()))
    println("  # length(workers())              : ", length(workers()))

    ## Set up logging directory
    logdir = set_logdir()
    println("- Logging directory: ", logdir)
    println()


    println("Experiment summarys:")
    println("  #problems:  ", length(problems))
    println("  algorithms: ", [a[:fnsolve_symbol] for a in algorithms])
    println("  seeds:      ", [a[:seeds] for a in algorithms])
    println()

    ## Logging object
    problem_to_algo = OrderedDict{String, Any}()
    problem_to_algo["problem_names"] = [pb[:pbname] for pb in problems]

    ## Run all algorithms once to precompile everything
    println("[", String(Dates.format(now(), "HHhMM SS")), "] Running algs once to precompile...")
    runallalgs()

    ## Solve
    for problem_descr in problems
        pb = problem_descr[:pb]
        pbname = problem_descr[:pbname]
        println("\n- [", String(Dates.format(now(), "HHhMM SS")), "] Dealing with new problem:")
        println(pb)
        @show problem_descr[:optimizer]
        @show problem_descr[:optimizer_params]

        
        algo_to_seedhist = OrderedDict{String, Any}()


        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        println("[", String(Dates.format(now(), "HHhMM SS")), "] Long global solve...")
        xsol = nothing
        fopt = nothing
        try
            # xsol = solve_progressivehedging(pb, ϵ_primal=1e-10, ϵ_dual=1e-10, maxtime=4*60*60, maxiter=1e6, printstep=10)
            xsol = problem_descr[:fnglobalsolve](pb)
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
            algo_to_seedhist[algo_descr[:algoname]][:seeds] = algo_descr[:seeds]
            merge!(algo_to_seedhist[algo_descr[:algoname]], algo_descr[:fnsolve_params])
            algo_to_seedhist[algo_descr[:algoname]][:optimizer] = problem_descr[:optimizer]
            algo_to_seedhist[algo_descr[:algoname]][:optimizer_params] = problem_descr[:optimizer_params]
            algo_to_seedhist[algo_descr[:algoname]][:nworkers] = length(workers())

            fnsolve = eval(algo_descr[:fnsolve_symbol])

            for seed in algo_descr[:seeds]
                println("  - [", String(Dates.format(now(), "HHhMM SS")), "] Solving for seed $seed")

                ## Set up log object
                hist = OrderedDict{Symbol, Any}()
                !isnothing(fopt) && (hist[:fopt] = fopt)
                !isnothing(xsol) && (hist[:approxsol] = xsol)

                println("\n--------------------------------------------------------")
                fnsolve(pb; algo_descr[:fnsolve_params]...,
                            optimizer = problem_descr[:optimizer],
                            optimizer_params = problem_descr[:optimizer_params],
                            hist = hist, 
                            seed = seed)
                println("--------------------------------------------------------\n")

                # Log seed, algo
                haskey(hist, :approxsol) && delete!(hist, :approxsol)
                algo_to_seedhist[algo_descr[:algoname]][seed] = hist
            end
        end

        problem_to_algo[pbname] = algo_to_seedhist
        problem_to_algo[pbname]["nstages"] = pb.nstages
        problem_to_algo[pbname]["nscenarios"] = pb.nscenarios
        problem_to_algo[pbname]["subpbdim"] = sum(length.(pb.stage_to_dim))
        # problem_to_algo[pbname]["nworkers"] = length(workers())
    end

    ## Write log information
    println("[", String(Dates.format(now(), "HHhMM SS")), "] All computations are completed. Writing logs.")
    writelogs(problem_to_algo, logdir)

    println("[", String(Dates.format(now(), "HHhMM SS")), "] All done.")
    return
end
