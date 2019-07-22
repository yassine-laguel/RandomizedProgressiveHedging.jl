# RPH.jl

*Randomized Progressive Hedging.*

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yassine-laguel.github.io/RPH.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yassine-laguel.github.io/RPH.jl/dev) | [![](https://api.travis-ci.com/yassine-laguel/RPH.jl.svg?token=aVdshbY2sAXsS8EzvkVr&branch=master)](https://travis-ci.com/yassine-laguel/RPH.jl) |

## Installation

```julia
]add https://github.com/yassine-laguel/RPH.jl#master
```

## Example

Many example scripts are available in the `example` folder. A good place to start is:
- the [doc](---)!
- the [`examples/tutorial.ipynb`](https://github.com/yassine-laguel/RPH.jl/blob/master/examples/tutorial.ipynb) jupyter notebook
- the script [`examples/simple_example.jl`](https://github.com/yassine-laguel/RPH.jl/blob/master/examples/tutorial.jl)

For distributed solve, launch julia as `julia -p 3` for 2 workers and a master thread.

## TODOs

### elementary
- [ ] [doc] Write documentation, readme
- [ ] [doc] Document functions properly
- [ ] [dev] make hydrothermal example fully parametric in nstages and nbraching
- [ ] Remove useless deps
- [ ] Refactor CVar