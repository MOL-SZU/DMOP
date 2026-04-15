# Scalable DMO-CNO Benchmark (PlatEMO-based)

This repository contains the code accompanying the paper:

**“A Scalable Benchmark Test Suite for Dynamic Multi-Objective Optimization with a Changing Number of Objectives.”**

The implementation is **based on PlatEMO** and focuses on **dynamic multi-objective optimization with a changing number of objectives (DMO-CNO)**.

## Background

Existing DMO-CNO benchmarks often replace a static objective count `m` with a time-varying `m(t)` directly inside DTLZ/WFG formulations. Because those formulations themselves depend on `m`, this unintentionally changes objective definitions over time.

This repository implements the paper’s key idea:

- Define a problem with a fixed maximum objective set.
- At each time step, activate only a subset of objectives.
- Keep objective functions themselves unchanged while only changing the active objective subset.

To avoid degeneracy issues from classical formulations, this benchmark uses **Minus-DTLZ** and **Minus-WFG** style objectives (all objectives mutually conflicting).

## What is implemented

### 1) Dynamic benchmark problem: `D_MIXED_SEQ.m`

`D_MIXED_SEQ` implements a stage-wise dynamic benchmark in which:

- a maximum objective pool is defined,
- active objective indices switch according to predefined stage sequences,
- environment changes occur at fixed frequencies,
- objective values are evaluated from DTLZ/WFG-style base constructions.

This realizes the paper’s “fixed objectives + dynamic active subset” benchmark design.

### 2) Performance metric: `MHV_Strict.m`

`MHV_Strict` computes a normalized hypervolume score for each stage:

- `HV_Pop`: hypervolume of algorithm population,
- `HV_PF`: hypervolume of reference Pareto front,
- `score = HV_Pop / HV_PF`.

This follows the paper’s motivation of making performance values comparable under different objective dimensions.

### 3) Experiment driver: `RunExperiments_WithMHV.m`

`RunExperiments_WithMHV` runs repeated experiments and records stage-wise results.

The included algorithms are:

- `DTAEA`
- `KTDMOEA`
- `LEC`
- `STA`

Algorithm implementations are under `Algorithms/` and are integrated in a PlatEMO-style workflow.

## Repository structure

- `D_MIXED_SEQ.m` – dynamic benchmark definition.
- `MHV_Strict.m` – MHV computation.
- `RunExperiments_WithMHV.m` – experiment script.
- `Algorithms/` – compared algorithms used in the paper.

## Usage

1. Prepare a MATLAB environment with the required **PlatEMO** base available.
2. Open MATLAB in this repository root.
3. Run:

```matlab
RunExperiments_WithMHV
```

Results are saved to `.mat` files as configured in the script.

## Notes

- This repo is an experimental research codebase built on PlatEMO conventions.
- Parameter settings (population size, frequency of change, number of runs, and objective-index sequences) should be configured in `RunExperiments_WithMHV.m` to reproduce specific paper settings.
