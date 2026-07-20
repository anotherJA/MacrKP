# code_MacrKP.f90

This is a Fortran 90 agent-based simulation I wrote to model the aggregation, growth, and fusion dynamics of a population of cells ("KP" cells) and macrophages ("Macr") confined in a spherical well. Each cell/macrophage (and any aggregate they form) is represented as a sphere with a position, radius, mass, and friction coefficient, evolving under overdamped (friction-dominated) dynamics with adhesion forces, gravity/confinement, and thermal-like noise.

It is used to explore how macrophage infiltration and adhesion affect the growth and packing of cell aggregates over time, with growth itself driven by a nutrient-diffusion-limited model (inside vs. outside diffusion coefficients, a critical nutrient concentration, and a growth timescale).

## What the code does

For a chosen initial number of cells and macrophages (`nc_cells`, `nc_macr`, picked from `ncelltable`/`nmacrtable` via the `num` index), the program:

1. **Initializes** positions (randomly placed near the center of the well), radii, masses, and friction/adhesion bookkeeping (`initialization`, `initial_pos`).
2. **Runs the main time loop** (`time_process`) for `tmax` steps, each step consisting of:
   - `boucle2` — for every aggregate: compute friction (`friction_function_new`), forces and updated velocity (`new_velocity_force`), and the new position (`new_position_force`); also handle growth (`growth`) and macrophage expulsion from over-full aggregates (`expulsion_macr`).
   - `boucle_fusion` — check pairs of aggregates for contact/adhesion and merge them (`fusion_function`), or split aggregates that exceed a critical size (`fragmentation`).
   - Periodically writes a snapshot (position, radius, composition, adhesion strengths) to the output `.dat` file.
3. **Finalizes** (`final_routine`) by printing summary statistics (min/max inter-aggregate distance, run duration) once the loop ends.

### Physical/numerical ingredients

- **Confinement**: aggregates are kept inside a spherical well of radius `radius_well` (`= xmax`); `rotate` and `position_function_new` handle geometric repositioning when the well boundary or another aggregate is hit.
- **Friction & forces**: `friction_function_new` computes anisotropic (normal/parallel) friction near the well surface, and `new_velocity_force`/`new_position_force` integrate the equations of motion (drag + adhesion + noise) with timestep `dt`.
- **Growth**: `growth` implements a diffusion-limited nutrient-uptake model (parameters `Din`, `Dext`, `fcrit`, `finfty`, `lpenetration`, `tau_growth`) that determines how many cells an aggregate gains or loses per step; `r3`/`rnbin` are small helper functions (cube-root radius mapping, binomial sampling).
- **Fusion/fragmentation**: `fusion_function` merges two aggregates that come into adhesive contact (probability set by `proba_cell`/`proba_macr`, decaying with `timeadh`/`timeadhm`/`tdelayc`/`tdelaym`); `fragmentation` splits aggregates above `ncritical`.
- **Macrophage packing**: `expulsion_macr` removes macrophages from an aggregate once the macrophage packing fraction `phimm` is exceeded.
- **RNG**: `ran2` is the classic *Numerical Recipes* long-period generator, seeded from the wall-clock time (`date_and_time`).

## Parameters

Almost all physical parameters are set in the block right after the variable declarations at the top of the program (radii, frictions, adhesion probabilities, diffusion constants, packing fractions `phi0`/`phi1`/`phimm`, timestep `dt=0.02`, `tmax=8250`, etc.). The outer `do 51 i51=num,num` loop is a leftover parameter-sweep pattern: as written it only runs the single index `num=9` (i.e. `nc_cells=2500`, `nc_macr=5000`), but changing `num` (or the loop bounds) selects a different point in `ncelltable`/`nmacrtable`.

## Output

Each run writes a single space-delimited `.dat` file named from three sweep-derived integer tags:

```
data<variable_name0><variable_name1><variable_name2>.dat
```

(macrophage count, `friction_macr*100`, `deltaG0*100`, each zero-padded to 5 digits). Every `1/dt` steps it appends a block starting with the current time, followed by one line per live aggregate:

```
x   y   z   radius   id   cell_fraction   n_cells   n_macr   adh_strength_cell   adh_strength_macr   radius_from_cell_count
```

## Compiling & running

Single-file, no external dependencies beyond a Fortran compiler (uses only intrinsic `date_and_time`; no MPI/OpenMP directives):

```bash
gfortran -O2 -o macrkp code_MacrKP.f90
./macrkp
```

Runtime scales with `tmax` × (`nc_cells`+`nc_macr`)², since several loops (e.g. `final_routine`'s distance check, `boucle_fusion`) are pairwise over all aggregates — for the default `nc_cells+nc_macr` ~ 7500 this is the main cost driver.
