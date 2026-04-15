# D_MIXED_SEQ

## Introduction

`D_MIXED_SEQ` is a dynamic many-objective optimization benchmark designed for evaluating evolutionary algorithms in changing environments. It supports stage-wise switching of active objective subsets, allowing different objectives to become active at different stages of the optimization process.

The benchmark integrates standard **DTLZ** and **WFG** formulations and also supports stage-dependent freezing of selected decision variables. This makes it a flexible testbed for studying the adaptability, robustness, and tracking ability of dynamic multi-objective and many-objective optimization algorithms. 

## Features

- Stage-wise switching of active objective subsets
- Mixed benchmark construction based on **DTLZ** and **WFG**
- Optional freezing of selected decision variables at different stages
- Built-in objective evaluation, reference front generation, and visualization
- MHV-based performance evaluation for dynamic optimization experiments 

## Code Structure

### 1. `D_MIXED_SEQ`

`D_MIXED_SEQ` defines the dynamic benchmark problem used in this project. It controls the active objective subsets at each stage, the mapping between objective indices and problem types, the frozen decision variables, and the objective evaluation process based on DTLZ/WFG formulations.

It also provides functions for dynamic stage switching, objective calculation, optimum generation through `GetOptimum`, and visualization through `DrawObj`. 

### 2. `MHV_Strict`

`MHV_Strict` is the performance indicator used in this project. It computes:
- `HV_Pop`: the hypervolume of the current population
- `HV_PF`: the hypervolume of the reference Pareto front
- `score`: the ratio `HV_Pop / HV_PF`, used as the final MHV value. 

The function first extracts objective values from the population and the reference Pareto front, normalizes them, and then computes hypervolume using exact calculation in 2D, exact calculation in small-scale 3D cases, and Monte Carlo estimation in higher dimensions. 

### 3. `RunExperiments_WithMHV`

`RunExperiments_WithMHV` is the main experiment script. It defines the population size, number of runs, and stage lengths, constructs benchmark settings, runs multiple algorithms, records the MHV value at the last generation of each stage, and saves the results into `.mat` files.

The current script evaluates the following algorithms:
- `DTAEA`
- `KTDMOEA`
- `LEC`
- `STA`

## Usage

Run the following script in MATLAB:

```matlab
RunExperiments_WithMHV

## Notes
This project is implemented based on the PlatEMO framework. Please make sure the required algorithms and utility functions are available in the MATLAB path before running the experiments. 
