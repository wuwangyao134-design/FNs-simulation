# FNs-simulation

This repository provides the MATLAB implementation used in our paper for simulating a fog computing network and evaluating the performance of resource scheduling using multi-objective optimization techniques.

## 🧩 Features

- Simulates fog node scheduling and task allocation
- Supports evaluation using multiple performance metrics:
  - Hypervolume (HV)
  - Generational Distance (GD)
  - Inverted Generational Distance (IGD)
- Includes example `.mat` datasets for experimental scenarios

## 📂 File Overview

| File | Description |
|------|-------------|
| `S_slot.m` | Main simulation script |
| `EvaluateParticle.m` | Objective function evaluator |
| `FindAllFronts.m`, `FindNonDominated.m` | Pareto front computation |
| `calculateHV.m`, `calculateGD.m`, `calculateIGD.m` | Metric calculation |
| `*.mat` | Sample datasets for different configurations |

## 🛠️ Requirements

- MATLAB R2024b or later
- No special toolboxes required (unless noted in paper)

## 📋 How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/wuwangyao134-design/FNs-simulation.git
