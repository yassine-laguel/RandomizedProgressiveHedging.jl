using DataStructures, GLPK, LinearAlgebra
# using Distributed
# @everywhere using RPH
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.status()

@everywhere using JuMP, RPH

using Dates, DelimitedFiles, JSON

GLOBAL_LOG_DIR = joinpath("/", "bettik", "PROJECTS", "pr-cvar")


function set_logdir()
    ## Setting working directory    
    date = String(Dates.format(now(), "yyyy_mm_dd-HHhMM"))
    logdir = joinpath(GLOBAL_LOG_DIR, "rphrun_$date")

    @show logdir
    @show isdir(logdir)
    !isdir(logdir) && mkpath(logdir)
    return logdir
end

function writelogs(problem_to_algo, logdir)
    logpath = joinpath(logdir, "full_logs.json")
    println("Writing logs to: ", logpath)
    
    open(logpath, "w") do w
        JSON.print(w, problem_to_algo)
    end
    return
end


function main()
    println("Entering main()")
    println("- Date is: ", String(Dates.format(now(), "yyyy_mm_dd-HHhMM")))
    println("  Path is: ", pwd())
    println("- ENV[\"OAR_NODEFILE\"]: ", get(ENV, "OAR_NODEFILE", ""))
    println("  Available cores:       ", vec(readdlm(get(ENV, "OAR_NODEFILE", ""), String)))

    ## Set up logging directory
    logdir = set_logdir()
    println("- Logging directory: ", logdir)
    println()


    ## Build problems to be solved
    problems = []

    push!(problems, (
        pb = ..,
    ))

    ## Build algorithms & params used for solve
    algorithms = []

    push!(algorithms, (
        fnsolve_symbol = ..,
        tlim = nothing,
        itmax = nothing,
    ))

    ## Set number of seeds to be tried
    seeds = 1:5


    ## Logging object
    problem_to_algo = OrderedDict{Symbol, OrderedDict}()

    println("#problems:  ", length(algorithms))
    println("algorithms: ", [a.fnsolve_symbol for a in algorithms])
    println("seeds:      ", seeds)
    println()

    ## Solve
    for problem_descr in problems
        pb = problem_descr.pb
        println("- [", String(Dates.format(now(), "HHhMM")), "] Dealing with new problem:")
        println(pb)

        ## First, solve pb up to reasonable precision. Get optimal objective, solution.
        xsol, fopt = 0, 0

        algo_to_seedhist = OrderedDict{Any, OrderedDict}()

        ## Run algorithms x seed
        for algo_descr in algorithms, seed in seeds
            println("- [", String(Dates.format(now(), "HHhMM")), "] Running algo:")
            println(algo_descr.fnsolve_symbol)
            println("seed:           ", seed)
            
            hist = OrderedDict{Symbol, Any}()
            
            fnsolve = eval(algo_descr.fnsolve_symbol)
            println("\n--------------------------------------------")
            fnsolve(pb; tlim=algo_descr.tlim, itmax=algo_descr.itmax, hist=hist, xsol=xsol, fopt=fopt, seed=seed)
            println("--------------------------------------------\n")

            # Log seed, algo
            algo_to_seedhist[algo_descr] = (
                seed = seed,
                hist = hist
            )
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
