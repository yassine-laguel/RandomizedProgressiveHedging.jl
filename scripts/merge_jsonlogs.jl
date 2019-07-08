using JSON, Gtk


function main()
    println("Merge jsons dealing with same problem, for different algorithms")
    
    # println("How many json should be merged ?")
    # njson = parse(Int, readline())

    # filepaths = String[]
    # for i in 1:njson
    #     push!(filepaths, open_dialog("Pick a json log"))
    # end

    filepaths = [ 
        "/home/gilles/Desktop/RPH.jl/logdir/RPH_num_exps/rphrun_2019_07_04-17h55/full_logs.json",
        "/home/gilles/Desktop/RPH.jl/logdir/RPH_num_exps/rphrun_2019_07_03-18h07/full_logs.json",
    ]
    @show filepaths

    logdicts = [JSON.parsefile(filepath) for filepath in filepaths]

    # println("-----")
    # for logdict in logdicts
    #     display(logdict)
    #     println("-----")
    # end


    ## Check same problem
    pbname = logdicts[1]["problem_names"][1]
    for logdict in logdicts
        @assert length(logdict["problem_names"]) == 1
        @assert logdict["problem_names"][1] == pbname
    end

    mergeddict = copy(logdicts[1])
    for logdict in logdicts[2:end]
        for (algname, algstat) in logdict[pbname]
            if !haskey(mergeddict[pbname], algname)
                merge!(mergeddict[pbname], Dict(algname=>algstat))
            else
                @warn "Not adding algorithm $algname, already present"
            end
        end
    end
    
    ## Write problem
    mergedjsonpath = save_dialog_native("Save merged json", GtkNullContainer(), String[])

    open(mergedjsonpath, "w") do w
        JSON.print(w, mergeddict)
    end
    println("All done.")
end

main()
