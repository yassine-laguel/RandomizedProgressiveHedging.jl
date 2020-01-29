# RandomizedProgressiveHedging.jl

*Randomized Progressive Hedging.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yassine-laguel.github.io/RandomizedProgressiveHedging.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yassine-laguel.github.io/RandomizedProgressiveHedging.jl/dev) | [![](https://api.travis-ci.com/yassine-laguel/RandomizedProgressiveHedging.jl.svg?token=aVdshbY2sAXsS8EzvkVr&branch=master)](https://travis-ci.com/yassine-laguel/RandomizedProgressiveHedging.jl) |

## Installation

The package si installed with the following command:
```julia
]add RandomizedProgressiveHedging.jl
```
GLPK, Ipopt are installed along as default solvers. Other solvers can of course be used, see [JuMP doc](http://www.juliaopt.org/JuMP.jl/v0.19.0/installation/#Getting-Solvers-1) for installation and `example/` scripts for use.

## Example

Many example scripts are available in the `example/` folder. A good place to start is:
- the [documentation](https://yassine-laguel.github.io/RandomizedProgressiveHedging.jl/stable)!
- the [`examples/tutorial.ipynb`](https://github.com/yassine-laguel/RandomizedProgressiveHedging.jl/blob/master/examples/tutorial.ipynb) jupyter notebook
- the script [`examples/simple_example.jl`](https://github.com/yassine-laguel/RandomizedProgressiveHedging.jl/blob/master/examples/tutorial.jl)

For distributed solve on local machine, launch e.g. julia as `julia -p 3` for 2 workers and a master thread. On a cluster, add workers with the package `Distributed` either `addprocs` on ssh connection, or use an adapted [cluster manager](https://github.com/JuliaParallel/ClusterManagers.jl).
