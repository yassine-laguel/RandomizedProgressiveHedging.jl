# RPH.jl

*Randomized Progressive Hedging.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yassine-laguel.github.io/RPH.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yassine-laguel.github.io/RPH.jl/dev) | [![](https://api.travis-ci.com/yassine-laguel/RPH.jl.svg?token=aVdshbY2sAXsS8EzvkVr&branch=master)](https://travis-ci.com/yassine-laguel/RPH.jl) |

## Installation

The package si installed with the following command:
```julia
]add https://github.com/yassine-laguel/RPH.jl#master
```
GLPK, Ipopt are installed along as default solvers. Other solvers can of course be used, see [JuMP doc](http://www.juliaopt.org/JuMP.jl/v0.19.0/installation/#Getting-Solvers-1) for installation and `example/` scripts for use.

## Example

Many example scripts are available in the `example/` folder. A good place to start is:
- the [doc](---)!
- the [`examples/tutorial.ipynb`](https://github.com/yassine-laguel/RPH.jl/blob/master/examples/tutorial.ipynb) jupyter notebook
- the script [`examples/simple_example.jl`](https://github.com/yassine-laguel/RPH.jl/blob/master/examples/tutorial.jl)

For distributed solve on local machine, launch e.g. julia as `julia -p 3` for 2 workers and a master thread. On a cluster, add workers with either `addprocs` on ssh connection, or use an adapted [cluster manager](https://github.com/JuliaParallel/ClusterManagers.jl).
