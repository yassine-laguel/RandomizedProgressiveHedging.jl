include("RPH.jl")

function makeproblem()
    n = 3

    ## Scenario 1
    model1 = Model()
    @variable model1 x[1:n]
    @NLobjective model1 Min sum((x[i]-1)^2 for i in 1:n)
    scenario1 = Scenario(model1)

    ## Scenario 2
    model2 = Model()
    @variable model2 x[1:n]
    @NLobjective model2 Min sum((x[i]-2)^2 for i in 1:n)
    scenario2 = Scenario(model2)

    ## Scenario 3
    model3 = Model()
    @variable model3 x[1:n]
    @NLobjective model3 Min sum((x[i]-3)^2 for i in 1:n)
    scenario3 = Scenario(model3)

    stageid_to_scenpart = [
        OrderedSet([BitSet(1:3)]),
        OrderedSet([BitSet(1), BitSet(2:3)]),
        OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),
    ]


    pb = Problem(
        [0.5, 0.25, 0.25],
        [scenario1, scenario2, scenario3],
        3,
        3,
        stageid_to_scenpart
    )

    return pb
end