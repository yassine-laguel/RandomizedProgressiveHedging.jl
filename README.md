# RPH.jl

Randomized Progressive Hedging.

## Installation

```julia
]add https://github.com/yassine-laguel/RPH.jl#master
]instantiate
]activate .
```

## Example

See scripts `exmaples/simple_example.jl`, `examples/hydrothermal_scheduling.jl` ~or notebook `hyforthermal_scheduling.ipynb`~. For distributed solve, launch julia as `julia -p 3` for 3 workers and a master thread.

## Questions

- What print info / log info ?
- Discuss stopping criteria
- Discuss algo parameters (upper delay bound, ...)

## TODOs

### Week 11->14 week
- [ ] Decide if the algorithms present now (vanilla, randomized, and asynchronous parallel) are the only ones we will consider. Typically, do/can we easily add a "distributed" version where the scenarios are local to the workers instead of sampled ?
- [ ] Create a more realistic test set (see the refs in the num. sec. of the paper) and proofread the code
- [ ] Start thinking about documentation, packaging, and creation of synthetic problems of different sizes

### Mid-term (june 26th)
- [ ] Run experiments and present nice figures/plots
- [ ] Package the code

### elementary
- [ ] [doc] Write documentation, readme
- [ ] [doc] Document functions properly
- [ ] [dev] make hydrothermal example fully parametric in nstages and nbraching
- [ ] [dev] make PH subpb solves distributed
- [ ] [dev] rph async - make sub pb solve parametric in solver
- [ ] [exp] write script/files for numerical experiments: first a direct/ph solve, then an rph solve with same time budget. Solution comparison, 'suboptimality' evolution, (smoothed ?) stepsize evolution. Max delay over time ?
- [ ] [?] Clean up code, stick to paper notations
- [ ] [?] Q sampling : parametrize & default to p

## Distributed computing with OAR

The 'SSHManager' doesnot work as is with the oar manager. Indeed, connection between avalaible nodes of a job are done with, e.g., `oarsh luke43`.

Two possibilities:
- write a `OarClusterManager`. Already considered by someone of imag, see [this discussion](https://discourse.julialang.org/t/which-workflow-to-launch-jobs-on-a-cluster/12532). `SSHManager` could be adapted, see [code](https://github.com/JuliaLang/julia/blob/55e36cc308b66d3472990a06b2797f9f9154ea0a/stdlib/Distributed/src/managers.jl#L5).
- `oarsh` is based on `ssh`, with specific parameters. Look into this to properly configure `SSHManager` - [doc here](https://www.grid5000.fr/w/Advanced_OAR#sharing_keys_between_jobs).