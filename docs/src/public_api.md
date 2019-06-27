# Public reference

## Building a problem
```@docs
RPH.ScenarioTree
RPH.Problem
```

## Algebra...
```@docs
objective_value
```

## Solving a problem
```@docs
solve_direct(pb::Problem; solver = with_optimizer(Ipopt.Optimizer))
solve_progressivehedging(pb::Problem)
solve_randomized_sync(pb::Problem)
solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario
```

## CVar extension
```@docs
cvar_problem(pb::Problem, cvar::CVar)
```