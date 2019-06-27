# Internal reference

## Scenario tree
```@docs
RPH.STreeNode

```

## Projections
```@docs
RPH.get_averagedtraj
RPH.get_averagedtraj!
RPH.nonanticipatory_projection
RPH.nonanticipatory_projection!
```

## Progressive hedging
```@docs
RPH.ph_subproblem_solve
```

## Synchronous
```@docs
RPH.randomizedsync_initialization!
RPH.randomizedsync_subpbsolve
```

## Asynchronous
```@docs
RPH.SubproblemTask
RPH.do_remote_work
RPH.randomizedasync_initialization!
RPH.init_workers
RPH.terminate_workers
```

## Others
```@docs
RPH.get_neighbydepth
RPH.get_scenariodim
RPH.dot
```