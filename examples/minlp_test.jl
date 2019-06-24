# using Juniper, Ipopt, JuMP

# using Cbc
# optimizer = Juniper.Optimizer
# params = Dict{Symbol,Any}()
# params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
# params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)

# using LinearAlgebra # for the dot product
# m = Model(with_optimizer(optimizer, params))

# v = [10,20,12,23,42]
# w = [12,45,12,22,21]
# @variable(m, x[1:5], Int)

# @objective(m, Max, dot(v,x))

# @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

# optimize!(m)

# # retrieve the objective value, corresponding x values and the status
# println(JuMP.value.(x))
# println(JuMP.objective_value(m))
# println(JuMP.termination_status(m))