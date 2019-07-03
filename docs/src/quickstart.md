# Quick start guide

This section aims provides an explanation of how to build and solve a problem using RPH.jl by solving a toy problem. The equivalent script and ijulia notebook can be found in the `example` folder.

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
& \sum_{t=1}^T C e_t + y_t \\
\text{s.t.}\quad
& q_t, y_t, e_t \ge 0 \\
& q_t \le W \\
& e_t+y_t \ge D \\
& q_1 = \bar{r}-y_1 \\
& q_t = q_{t-1}+r[\xi_t]-y_t, \; t = 2, \ldots, T.
\end{aligned}
```

where ``C = 5``, ``W = 8``, ``D = 6``, ``r = [2, 10]``. A scenario is defined by ``(\xi_t)_{t=2, \ldots, T}``, for ``\xi_t\in\{1,2\}``.

### Representing a scenario

A scenario is represented by the following structure:
```julia
struct HydroThermalScenario <: AbstractScenario
    weather::Vector{Int}
end
```
Here, the attribut `weather` will hold one realisation of ``(\xi_t)_{t=2, \ldots, T}``.

Along with this scenario structure, the function laying out the scenario objective function ``f_s`` needs to be defined.
It takes as input the JuMP model that will hold ``f_s``, an instance of the previously defined scenario structure, and the identifier of the scenario.
```julia
function build_fs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    C = 5
    W = 8
    D = 6
    rain = [2, 10]

    T = length(s.weather)+1
    Y = @variable(model, [1:3*T], base_name="y_s$id_scen")

    q = [Y[1+3*k] for k in 0:T-1]
    y = [Y[2+3*k] for k in 0:T-1]
    e = [Y[3+3*k] for k in 0:T-1]

    ## State variables constraints
    @constraint(model, Y[:] .>= 0)      # positivity constraint
    @constraint(model, q .<= W)         # reservoir max capacity
    @constraint(model, e .+ y .>= D)    # meet demand
    
    ## Dynamic constraints
    @constraint(model, q[1] == sum(rain)/length(rain) - y[1])
    @constraint(model, [t=2:T], q[t] == q[t-1] - y[t] + rain[s.weather[t-1]+1])
    
    objexpr = C*sum(e) + sum(y)

    return Y, objexpr, []
end
```
!!! note

    - The last item returned by the function should be the reference of constraints used to build the objective, none here. Such constraints can appear when modelling a ``\max(u, v)`` in the objective as a variable ``m``, under the linear constraints ``m\ge u`` and ``m\ge v``.


### Representing the scenario tree
The scenario tree represents the stage up to which scenarios are equal.

If the problem scenario tree is a [perfect *m*-ary tree](https://en.wikipedia.org/wiki/M-ary_tree#Types_of_m-ary_trees), it can be built using a buit-in function:
```julia
scenariotree = ScenarioTree(; depth=T, nbranching=2)
```

If the tree is not regular, or quite simple, it can be built by writing specifically the partition of equivalent scenarios per stage. A simple exmaple would be:
```julia
stageid_to_scenpart = [
    OrderedSet([BitSet(1:3)]),                      # Stage 1
    OrderedSet([BitSet(1), BitSet(2:3)]),           # Stage 2
    OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),  # Stage 3
]
```
!!! note
    However this method is not efficient, and whenever possible, builtin functions should be priviledged.

### Building the `Problem`


```julia
scenid_to_weather(scen_id, T) = return [mod(floor(Int, scen_id / 2^i), 2) for i in T-1:-1:0]

T = 5
nbranching = 2

p = 0.5

nscenarios = 2^(T-1)
scenarios = HydroThermalScenario[ HydroThermalScenario( scenid_to_weather(scen_id, T-1) ) for scen_id in 0:nscenarios-1]
probas = [ prod(v*p + (1-v)*(1-p) for v in scenid_to_weather(scen_id, T-1)) for scen_id in 1:nscenarios ]

stage_to_dim = [3*i-2:3*i for i in 1:T]
scenariotree = ScenarioTree(; depth=T, nbranching=2)


pb = Problem(
    scenarios::Vector{HydroThermalScenario},
    build_fs!::Function,
    probas::Vector{Float64},
    nscenarios::Int,
    T::Int,
    stage_to_dim::Vector{UnitRange{Int}},
    scenariotree::ScenarioTree,
)
```

## Solving a problem

### Explicitly solving the problem
```julia
y_direct = solve_direct(pb)
println("\nDirect solve output is:")
display(y_direct)
@show objective_value(pb, y_direct)
```

### Solving with Progressive Hedging
```julia
y_PH = solve_progressivehedging(pb, ϵ_primal=1e-4, ϵ_dual=1e-4, printstep=5)
println("\nSequential solve output is:")
display(y_PH)
@show objective_value(pb, y_PH)
```

### Solving with Randomized Progressive Hedging
```julia
y_sync = solve_randomized_sync(pb, maxtime=5, printstep=50)
println("\nSynchronous solve output is:")
display(y_sync)
@show objective_value(pb, y_sync)
```

### Solving with Parallel Randomized Progressive Hedging
```julia
y_par = solve_randomized_par(pb, maxtime=5, printstep=50)
println("\nSynchronous solve output is:")
display(y_par)
@show objective_value(pb, y_par)
```

### Solving with Asynchronous Randomized Progressive Hedging
Asynchronous solve leverages the distributed capacities of julia. In order to be used, workers need to be available. Local or remote workers can be added with [`addprocs`](https://docs.julialang.org/en/v1/stdlib/Distributed/#Distributed.addprocs).

`RPH` and `JuMP` packages need to be available for all workers, along with the scenario object and objective build function.

```julia
addprocs(3)     # add 3 workers besides master
@everywhere using RPH, JuMP

@everywhere HydroThermalScenario, build_fs!

y_async = solve_randomized_async(pb, maxtime=5, printstep=100)
println("Asynchronous solve output is:")
display(y_async)
@show objective_value(pb, y_par)
```
