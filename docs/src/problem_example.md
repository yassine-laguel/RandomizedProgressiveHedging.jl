# Tutorial

Detail HydroWater scheduling example here.

## Building the problem

The problem is:

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

## Solving with Progressive Hedging

## Solving with Randomized Progressive Hedging
