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

        @printf " - %-16s  %-6i\n" "maxiter" maxiter
        @printf " - %-16s  %-6i\n" "maxtime (s)" maxtime
        for (param, pval) in kwargs
            @printf " - %-16s  " param
            println(pval)
        end
    end
    return
end


function rand_getresidual(pb, z, z_old)
    return norm(pb, z-z_old)
end

function rand_hasconverged(pb, z, residual, ϵ_abs, ϵ_rel)
    return residual < ϵ_abs + ϵ_rel * norm(pb, z)
end
