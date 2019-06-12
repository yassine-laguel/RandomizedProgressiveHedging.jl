# RPH.jl Documentation

Hello there. Presentation of the module

```@contents
Pages = ["index.md", "problem_example.md"]
```

## Building a problem
lorem...

## Solving a problem
```@docs
solve_direct(pb::Problem; solver = with_optimizer(Ipopt.Optimizer))
solve_progressivehedging(pb::Problem)
solve_randomized_sync(pb::Problem)
solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario
```
