# &#x20;CO2 Diffusion Model

This repository contains a cleaned MATLAB package for reproducing the local concentration results from the diffusion-reaction model.

## Files

* `run\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_diffusion\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_model.m`  
Main diffusion model program.
* `README.md`  
Run instructions and repository usage.

## What is included

* Main model algorithm for local concentration simulation
* PDE-based diffusion-reaction solver
* Clean parameter-handling workflow

## Required input files

Place your own cleaned input files in a `data/` folder:

* `data/bulk.xlsx`
* `data/partial\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_current.xlsx`

### Expected format

`bulk.xlsx`

* Column 1: pressure
* Column 3: bulk CO2 concentration
* Column 5: bulk HCO3- concentration
* Column 6: bulk CO3^2- concentration
* Column 8: bulk OH- concentration

`partial\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_current.xlsx`

* Each sheet corresponds to one voltage / operating condition.
* Column 1: voltage
* Column 2: HAC current density
* Column 3: H2 current density

## How to run

```matlab
results = run\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_diffusion\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_model();
```

The script reads input files from the `data/` folder and saves results to `results/diffusion\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_model\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_results.mat`.

## Notes

* The model uses `pdepe`.
* The diffusion layer thickness and kinetic constants are defined inside the configuration block.
* For a GitHub release, keep only the code files and a short README.

