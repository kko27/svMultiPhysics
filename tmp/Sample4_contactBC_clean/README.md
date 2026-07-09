This is a cleaned draft of the original `Sample4` rigid-plane contact setup.

Key changes from the original case:
- Rigid-plane penalty contact is applied only on `top_cap`.
- `bottom_cap` is explicitly constrained in all three directions using
  `Effective_direction (1, 1, 1)` to remove rigid-body modes.
- The displacement history file is replaced with a solver-aligned scalar
  time history over the actual simulation interval `[0, 10]`.
- The nonlinear and linear solver settings are made closer to the stable
  in-repo rigid-plane contact example.
