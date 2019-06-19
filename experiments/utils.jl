
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

    return
end


function runallalgs()
    pb = build_simpleexample()

    y_direct = solve_direct(pb)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    y_PH = solve_progressivehedging(pb)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    y_synch = solve_randomized_sync(pb)
    println("\nSynchronous solve output is:")
    display(y_synch)
    
    y_asynch = solve_randomized_async(pb)
    println("Asynchronous solve output is:")
    display(y_asynch)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_synch)
    @show norm(y_direct - y_asynch)

    return
end
