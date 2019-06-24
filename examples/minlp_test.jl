using Juniper, Ipopt, JuMP
using Cbc
using LinearAlgebra # for the dot product

function main()
    optimizer = Juniper.Optimizer
    params = Dict{Symbol,Any}()
    params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
    params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)

    model = Model(with_optimizer(optimizer, params))

    id_scen = 2
    B = 3
    T = 2

    c_H = [1 for b in 1:B]       # dams electiricity prod costs
    c_E = 6             # externel elec buy cost
    D = 6               # demand at each time step
    W = 8               # dam water capacity
    W1 = 6              # initial state of dams
    rain = [2, 10]

    # Convert weathertype::Int into stage to rain level
    stage_to_rainlevel = [1, 2]

    qs = @variable(model, [1:B*T], base_name="qs_s$id_scen")
    ys = @variable(model, [1:B*T], base_name="ys_s$id_scen", Int)
    e = @variable(model, [1:T], base_name="e_s$id_scen")

    # positivity constraint
    for t in 1:T
        @constraint(model, e[t] >= 0)
    end
    @constraint(model, qs .>= 0)
    @constraint(model, ys .>= 0)


    # Meet demand
    @constraint(model, [t=1:T], sum(ys[(t-1)*B+1:t*B]) + e[t] >= D)

    # Reservoir max capacity
    @constraint(model, qs .<= W)

    ## Dynamic constraints
    @constraint(model, qs[1:B] .== W1 - ys[1:B])
    @constraint(model, [t=2:T], qs[(t-1)*B + 1:(t-1)*B + B] .== qs[(t-2)*B + 1:(t-2)*B + B] - ys[(t-1)*B + 1:(t-1)*B + B] .+ rain[stage_to_rainlevel[t]])

    # objexpr = @NLexpression(model, sum(sum(c_H[i]*ys[(t-1)B + i] for i in 1:B) + c_E * e[t] for t in 1:T))
    # @NLobjective model Min objexpr

    objexpr = @expression(model, sum(sum(c_H[i]*ys[(t-1)B + i] for i in 1:B) + c_E * e[t] for t in 1:T) + ys[1]^2)
    @objective model Min objexpr

    println(model)


    @show model.nlp_data
    @show typeof(model)
    #, EmptyNLPEvaluator

    # Y = collect(Iterators.flatten([ union(ys[(t-1)*B+1:t*B], qs[(t-1)*B+1:t*B], e[t]) for t in 1:T] ))

    # v = [10,20,12,23,42]
    # w = [12,45,12,22,21]
    # @variable(m, x[1:5], Int)

    # @objective(m, Max, dot(v,x))

    # @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    optimize!(m)

    # retrieve the objective value, corresponding x values and the status
    println(JuMP.value.(x))
    println(JuMP.objective_value(m))
    println(JuMP.termination_status(m))
end

main()