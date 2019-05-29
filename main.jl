# using RPH
using Ipopt
include("src/RPH.jl")
include("src/testcases.jl")

function main()
    pb = makeproblem()
    
    model1 = Model(with_optimizer(Ipopt.Optimizer))
    @variable model1 x
    model1_obj = @NLexpression model1 (x-2)^2 + 1/3
    @NLobjective model1 Min (x-2)^2 + 1/3

    optimize!(model1)
    @show value(x)

    
    model = Model(with_optimizer(Ipopt.Optimizer))


    return pb
end
main()
