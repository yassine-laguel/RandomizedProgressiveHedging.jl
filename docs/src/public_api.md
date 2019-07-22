# Public reference

## Building a problem
```@docs
AbstractScenario
build_fs!
RPH.ScenarioTree
RPH.Problem
```

## Solving a problem
```@docs
solve_direct(pb::Problem; solver = with_optimizer(Ipopt.Optimizer))
solve_progressivehedging(pb::Problem)
solve_randomized_sync(pb::Problem)
solve_randomized_par
solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario
```

## Other
```@docs
objective_value
cvar_problem
```