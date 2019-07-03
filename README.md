# RPH.jl

*Randomized Progressive Hedging.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yassine-laguel.github.io/RPH.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yassine-laguel.github.io/RPH.jl/dev) | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

[![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url]

## Installation

```julia
]add add https://github.com/GillesBareilles/OarClusterManager.jl#master
]add https://github.com/yassine-laguel/RPH.jl#master
]instantiate
]activate .
```

## Example

See scripts `exmaples/simple_example.jl`, `examples/hydrothermal_scheduling.jl` ~or notebook `hyforthermal_scheduling.ipynb`~. For distributed solve, launch julia as `julia -p 3` for 3 workers and a master thread.

## TODOs

### Week 11->14 week
- [ ] Create a more realistic test set (see the refs in the num. sec. of the paper) and proofread the code

### Mid-term (june 26th)
- [ ] Run experiments and present nice figures/plots
- [ ] Package the code

### elementary
- [ ] [doc] Write documentation, readme
- [ ] [doc] Document functions properly
- [ ] [dev] make hydrothermal example fully parametric in nstages and nbraching

## Distributed computing with OAR

The 'SSHManager' doesnot work as is with the oar manager. Indeed, connection between avalaible nodes of a job are done with, e.g., `oarsh luke43`.

Use `OarClusterManager`, at (https://github.com/GillesBareilles/OarClusterManager.jl#master).