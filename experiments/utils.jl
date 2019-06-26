
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

    # ## Copy stdoud, stderr file redirect of oar to right folder
    # if haskey(ENV, "OAR_JOB_ID")
    #     for ext in ["stdout", "stderr"]
    #         filename = filter(x->occursin("$(ENV["OAR_JOB_ID"]).$ext", x), readdir(homedir()))
    #         if length(filename) == 1
    #             cp(joinpath(homedir(), filename[1]), joinpath(logdir, filename[1]))
    #         else
    #             println("No $ext file found.")
    #         end
    #     end
    # end

    return
end


function runallalgs(; runasync=true)
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
    
    runasync && print("solve_randomized_async... ")
    runasync && (y_asynch = solve_randomized_async(pb; printlev=0, maxtime=0.5))
    runasync && println("done.")
    # println("Asynchronous solve output is:")
    # display(y_asynch)

    # @show norm(y_direct - y_PH)
    # @show norm(y_direct - y_synch)
    # @show norm(y_direct - y_asynch)

    return
end
