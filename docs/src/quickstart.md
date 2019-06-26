# Quick start guide

## Installation
RPH.jl is a pure julia package. It can be installed from julia by using the built-in package manager:
```julia
using Pkg
Pkg.add("https://github.com/yassine-laguel/RPH.jl")
```

## Getting solvers
RPH depends on other solvers to optimize the subproblems. All solvers interfaced with JuMP, the julia mathematical programming language, can be used in RPH, a list of which can be found at [JuMP's documentation](http://www.juliaopt.org/JuMP.jl/v0.19.0/installation/#Getting-Solvers-1).

Note that all algorithms layout subproblem with quadratic objectives. Default subproblem solver is the interior point algorithm Ipopt.

## Laying out a problem
We take the following problem as example:

```math
\begin{aligned}
\underset{x}{\text{minimize}}\quad
& \sum_{t=1}^T C p_t + y_t \\
\text{s.t.}\quad
& x_t, y_t, p_t \ge 0 \\
& x_t \le W \\
& p_t+y_t \ge D \\
& x_1 = \bar{r}-y_1 \\
& x_t = x_{t-1}+r[\xi_t]-y_t, \; t = 2, \ldots, T.
\end{aligned}
```

### Representing a scenario

A scenario is represented by the following structure:
```julia
struct HydroThermalScenario <: RPH.AbstractScenario
    weather::Vector{Int}
end
```

Along with this scenario structure, the function laying out the scenario objective function ``f_s`` needs to be defined.
It takes as input the JuMP model that will hold ``f_s``, an instance of the previously defined scenario structure, and the identifier of the scenario.
```julia
function build_fs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    C = 5
    W = 8
    D = 6
    rain = [2, 10]

    nstages = length(s.weather)+1
    Y = @variable(model, [1:3*nstages], base_name="y_s$id_scen")

    x = [Y[1+3*k] for k in 0:nstages-1]
    y = [Y[2+3*k] for k in 0:nstages-1]
    p = [Y[3+3*k] for k in 0:nstages-1]

    ## State variables constraints
    @constraint(model, Y[:] .>= 0)      # positivity constraint
    @constraint(model, x .<= W)         # reservoir max capacity
    @constraint(model, p .+ y .>= D)    # meet demand
    
    ## Dynamic constraints
    @constraint(model, x[1] == sum(rain)/length(rain) - y[1])
    @constraint(model, [t=2:n], x[t] == x[t-1] - y[t] + rain[stage_to_rainlevel[t]])
    
    objexpr = C*sum(p) + sum(y)

    return Y, objexpr, []
end
```
!!! note

    - The last item returned by the function should be the reference of constraints used to build the objective, none here. Such constraints can appear when modelling a ``\max(u, v)`` in the objective as a variable ``m``, under the linear constraints ``m\ge u`` and ``m\ge v``.


### Representing the scenario tree
The scenario tree represents the stage up to which scenarios are equal.

It can be built by writing specifically the partition of scenarios per stage. A simple exmaple would be:
```julia
stageid_to_scenpart = [
    OrderedSet([BitSet(1:3)]),                      # Stage 1
    OrderedSet([BitSet(1), BitSet(2:3)]),           # Stage 2
    OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),  # Stage 3
]
```

However this method is not efficient, and when possible, builtin functions should be priviledged. When the scenario tree has a known depth and each node has the same number of childs, one should prefer:
```julia
scenariotree = ScenarioTree(; depth=nstages, nbranching=2)
```

### Building the `Problem`

```julia
Problem(
    scenarios,
    build_fs!,
    probas,
    nscenarios, 
    nstages,
    dim_to_subspace,
    scenariotree
)
```

## Solving a problem

### Solving with Progressive Hedging

### Solving with Randomized Progressive Hedging

### Solving with Asynchronous Randomized Progressive Hedging
