
function set_logdir()
    ## Setting working directory    
    date = String(Dates.format(now(), "yyyy_mm_dd-HHhMM"))
    logdir = joinpath(GLOBAL_LOG_DIR, "rphrun_$date")

    !isdir(logdir) && mkpath(logdir)
    return logdir
end

function writelogs(problem_to_algo, logdir)
    logpath = joinpath(logdir, "full_logs.json")
    println("Writing logs to: ", logpath)
    
    open(logpath, "w") do w
        JSON.print(w, problem_to_algo)
    end

    ## Copy stdoud, stderr file redirect of oar to right folder
    for ext in ["stdout", "stderr"]
        filename = filter(x->occursin("20102139.$ext", x), readdir())
        if length(filename) == 1
            cp(filename[1], joinpath(logdir, filename[1]))
        end
    end

    isfile()
    return
end

function checknbcores(algorithms)
    nbavailcore = length(vec(readdlm(get(ENV, "OAR_NODEFILE", ""), String)))
    for a in algorithms
        !isnothing(a.ncores) && @assert a.ncores <= nbavailcore "Algorithm $(a.algoname) requires $(a.ncores), $nbavailcore available."
    end
    return
end

function allocateworkers(ncores)
    @show ncores, workers(), Vector([1])
    workers() !== Vector([1]) && (rmprocs(workers()); println("removing workers"))
    if !isnothing(ncores)
        ncores_masternode = get_ncoresmaster()
        remotenodecores = get_remotehosts()

	@show workers()
        addprocs(ncores_masternode-1)
	@show workers()
	@show remotenodecores
	@show remotenodecores[1:ncores-ncores_masternode]
        addprocs_oar(remotenodecores[1:ncores-ncores_masternode])
	@show workers()
        sleep(0.1)
        @assert length(workers()) == ncores-1

	## ensure package availability...
	# @everywhere eval(Expr(:using, :Pkg))
	# @everywhere Pkg.activate(".")
	@everywhere push!(LOAD_PATH, joinpath(homedir(), "RPH.jl"))
	@everywhere eval(Expr(:using, :RPH))
	@everywhere eval(Expr(:using, :JuMP))
    end
    return
end

function runallalgs()
    allocateworkers(length(vec(readdlm(get(ENV, "OAR_NODEFILE", ""), String))))

    pb = build_simpleexample()

    # y_direct = solve_direct(pb)
    println("\nDirect solve output is:")
    # display(y_direct)
    println("")

    # y_PH = solve_progressivehedging(pb)
    println("\nSequential solve output is:")
    # display(y_PH)
    println("")

    # y_synch = solve_randomized_sync(pb)
    println("\nSynchronous solve output is:")
    # display(y_synch)
    
    y_asynch = solve_randomized_async(pb)
    println("Asynchronous solve output is:")
    display(y_asynch)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_synch)
    @show norm(y_direct - y_asynch)

    return
end
