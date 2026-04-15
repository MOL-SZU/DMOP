# D_MIXED_SEQ

`D_MIXED_SEQ` is a dynamic many-objective optimization benchmark for evaluating evolutionary algorithms in changing environments. It supports stage-wise objective switching, where different subsets of objectives become active over time. :contentReference[oaicite:0]{index=0}

The benchmark combines standard **DTLZ** and **WFG** formulations and also allows selected decision variables to be frozen at different stages. This provides a flexible testbed for studying the adaptability and robustness of dynamic multi-objective and many-objective optimization algorithms. :contentReference[oaicite:1]{index=1}

## Features

- Dynamic switching of active objective subsets across optimization stages
- Support for mixed benchmark construction using **DTLZ** and **WFG**
- Stage-dependent freezing of selected decision variables
- Built-in objective evaluation, optimum generation, and visualization
- Suitable for dynamic multi-objective and many-objective optimization studies :contentReference[oaicite:2]{index=2}

## Configuration

The benchmark can be configured through the following parameters:

- `taut`: number of generations per stage
- `SeqString`: sequence of active objective subsets
- `ObjConfig`: mapping between objective indices and benchmark types
- `FrozenStr`: stage-wise specification of frozen decision variables :contentReference[oaicite:3]{index=3}

## Usage

`D_MIXED_SEQ` is implemented as a problem class and can be used within the optimization framework as a dynamic test problem. During evolution, the benchmark updates the active objective set according to the current stage and evaluates solutions based on the corresponding DTLZ or WFG formulation. :contentReference[oaicite:4]{index=4}

## Purpose

This benchmark is designed to help researchers study how optimization algorithms respond to nonstationary environments, especially in terms of convergence, diversity maintenance, robustness, and tracking ability. :contentReference[oaicite:5]{index=5}
