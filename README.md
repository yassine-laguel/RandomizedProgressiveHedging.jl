# RandomizedProgressiveHedging.jl

This module aims to solve multistage stochastic problems by randomized versions of the progressive hedging algorithm. It comes with the [companion paper](https://hal.archives-ouvertes.fr/hal-02946615/document) published Annals of Operations Research.

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yassine-laguel.github.io/RandomizedProgressiveHedging.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yassine-laguel.github.io/RandomizedProgressiveHedging.jl/dev) | [![](https://api.travis-ci.com/yassine-laguel/RandomizedProgressiveHedging.jl.svg?token=aVdshbY2sAXsS8EzvkVr&branch=master)](https://travis-ci.com/yassine-laguel/RandomizedProgressiveHedging.jl) |

## Installation

The package is installed with the following command:
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

## Authors
[Gilles Bareilles](https://gbareilles.fr)  
[Yassine Laguel](https://yassine-laguel.github.io)  
[Dmitry Grischenko](https://grishchenko.org)  
[Franck Iutzeler](http://www.iutzeler.org)  
[Jérôme Malick](https://ljk.imag.fr/membres/Jerome.Malick/)  

## Cite

If you found this package useful, please cite the following work.

```
@article{bareilles:hal-02946615,
  TITLE = {{Randomized Progressive Hedging methods for Multi-stage Stochastic Programming}},
  AUTHOR = {Bareilles, Gilles and Laguel, Yassine and Grishchenko, Dmitry and Iutzeler, Franck and Malick, J{\'e}r{\^o}me},
  URL = {https://hal.archives-ouvertes.fr/hal-02946615},
  JOURNAL = {{Annals of Operations Research}},
  PUBLISHER = {{Springer Verlag}},
  YEAR = {2020},
  PDF = {https://hal.archives-ouvertes.fr/hal-02946615/file/main.pdf},
  HAL_ID = {hal-02946615},
  HAL_VERSION = {v1},
}
```
