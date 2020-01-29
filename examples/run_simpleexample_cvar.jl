using Distributed
@everywhere using RandomizedProgressiveHedging, JuMP

include("build_simpleexample.jl")
using GLPK, Ipopt

function main()
    pb = build_simpleexample()

    hist=OrderedDict{Symbol, Any}(
        :approxsol =>   [0.468975  1.72212  1.0      1.0
                         0.468975  1.72212  2.62547  2.0
                         0.468975  1.72212  2.62547  3.0]
    )

    println("Full problem is:")
    println(pb)

    cvar_lev = 0.8
    pbcvar = cvar_problem(pb, cvar_lev)


    function callback(cvar_pb::Problem, x, hist)
        @assert hist !== nothing

        !haskey(hist, :cvarvalues) && (hist[:cvarvalues]=Float64[])

        fvalues = [objective_value(pb, x[:, 2:end], id_scen) for id_scen in 1:pb.nscenarios]
        model = Model(with_optimizer(GLPK.Optimizer))

        @variable(model, eta)
        @variable(model, m[1:pb.nscenarios])
        @objective(model, Min, eta + 1/(1-cvar_lev) * sum(pb.probas[i] * m[i] for i in 1:pb.nscenarios))
        @constraint(model, m .>= 0 )
        @constraint(model, [i in 1:pb.nscenarios], m[i]>= fvalues[i] - eta)

        optimize!(model)

        eta_opt = JuMP.value(eta)
        push!(hist[:cvarvalues], eta_opt)
        @show(eta_opt)
    end

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pbcvar, optimizer = Ipopt.Optimizer, printlev=0)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    hist = OrderedDict{Symbol, Any}()
    y_PH = solve_progressivehedging(pbcvar, maxiter=150, maxtime=40, printstep=5, hist=hist, callback=callback)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    hist = OrderedDict{Symbol, Any}()
    y_sync = solve_randomized_sync(pbcvar, maxiter=600, maxtime=40, printstep=100, hist=hist, callback=callback)
    println("\nSynchronous solve output is:")
    display(y_sync)

    # #########################################################
    # ## Problem solve: asynchronous (parallelized) version of PH
    hist = OrderedDict{Symbol, Any}()
    y_async = solve_randomized_async(pbcvar, maxiter=1e5, maxtime=20, printstep=100, hist=hist, callback=callback)
    println("Asynchronous solve output is:")
    display(y_async)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_sync)
    @show norm(y_direct - y_async)

    return
end

main()
