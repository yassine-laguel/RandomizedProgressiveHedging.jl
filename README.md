# RPH.jl

Randomized Progressive Hedging.

## Installation

```julia
]add https://github.com/yassine-laguel/RPH.jl#dev_gilles
]instantiate
]activate .
```

## Example

See scripts `exmaples/simple_example.jl`, `examples/hydrothermal_scheduling.jl` ~or notebook `hyforthermal_scheduling.ipynb`~. For distributed solve, launch julia as `julia -p 3` for 3 workers and a master thread.

## Questions

- What initialization for x: 0, subpbs solutions ?
- Interfaces with JuMP world packages : StructJuMP ?
- What print info / log info ?

## TODOs

### Week 11->14 week
- [ ] Decide if the algorithms present now (vanilla, randomized, and asynchronous parallel) are the only ones we will consider. Typically, do/can we easily add a "distributed" version where the scenarios are local to the workers instead of sampled ?
- [ ] Create a more realistic test set (see the refs in the num. sec. of the paper) and proofread the code
- [ ] Start thinking about documentation, packaging, and creation of synthetic problems of different sizes

### Mid-term (june 26th)
- [ ] Run experiments and present nice figures/plots
- [ ] Package the code

### Random
- [ ] For fair comparison with RPH, make PH subpb computation also distributed.
- [ ] Clean up code, stick to paper notations

