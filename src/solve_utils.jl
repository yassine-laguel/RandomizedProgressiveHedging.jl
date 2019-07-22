function display_algopb_stats(pb::Problem, algoname, printlev; maxiter, maxtime, kwargs...)
    if printlev > 0
        println("\n--------------------------------------------------------")
        println("--- $algoname")
        println("--------------------------------------------------------")
        println("Problem with:")
        println(" - nb scenarios       : ", pb.nscenarios)
        println(" - nb stages          : ", pb.nstages)
        println(" - variable dimension : ", get_scenariodim(pb))
        println("Algorithm parameters:")

        @printf " - %-10s  %-6i    - %-10s  %-6i" "maxiter" maxiter "maxtime" maxtime
        for (param, pval) in kwargs
            @printf "\n - %-10s  " param
            display(pval)
        end
        println("")
    end
    return
end