
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

    return
end


function runallalgs()
    runasync=false
    workers() != Vector([1]) && (runasync=true)

    pb = build_simpleexample()

    # y_direct = solve_direct(pb; printlev=0, optimizer=Ipopt.Optimizer)
    # # println("\nDirect solve output is:")
    # # display(y_direct)
    # # println("")

    print("solve_progressivehedging... ")
    y_PH = solve_progressivehedging(pb; printlev=0, maxtime=0.5)
    println("done.")
    # println("\nSequential solve output is:")
    # display(y_PH)
    # println("")

    print("solve_randomized_sync... ")
    y_synch = solve_randomized_sync(pb; printlev=0, maxtime=0.5)
    println("done.")
    # println("\nSynchronous solve output is:")
    # display(y_synch)

    runasync && print("solve_randomized_par... ")
    runasync && (y_asynch = solve_randomized_par(pb; printlev=0, maxtime=20))
    runasync && println("done.")
    # println("Asynchronous solve output is:")
    # display(y_asynch)
    
    runasync && print("solve_randomized_async... ")
    runasync && (y_asynch = solve_randomized_async(pb; printlev=0, maxtime=20))
    runasync && println("done.")
    # println("Asynchronous solve output is:")
    # display(y_asynch)

    # @show norm(y_direct - y_PH)
    # @show norm(y_direct - y_synch)
    # @show norm(y_direct - y_asynch)

    return
end
