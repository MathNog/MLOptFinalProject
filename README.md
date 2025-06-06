# Extending Multistage Day-Ahead Scheduling with Piecewise Linear Decision Rules

## Author  
Matheus Nogueira, Matheus Alves and Luis Fernando Duarte

## About this Repository

This repository contains an academic implementation and methodological extension of the paper:

> Rodrigues, M., Street, A., & Arroyo, J.M. (2024).  
> *Multistage Day-Ahead Scheduling of Energy and Reserves*.  
> Electric Power Systems Research, 235, 110793.  
> [DOI: 10.1016/j.epsr.2024.110793](https://doi.org/10.1016/j.epsr.2024.110793)

The original paper proposes a **multistage stochastic programming model** for joint energy and reserve scheduling under uncertainty using **regularized Linear Decision Rules (LDR)**. In this project, we explore the possibility of **replacing the LDR by a Piecewise Linear Decision Rule (PLDR)** to improve modeling flexibility in the presence of nonlinear dynamics and structural changes in uncertainty.

### Objectives of the Extension

- Investigate the tractability and computational performance of implementing PLDR in place of traditional LDR.
- Evaluate the out-of-sample performance of PLDR-based policies using the same simulation framework proposed by the authors.
- Assess whether the increased representational power of PLDR justifies its computational cost in practical scenarios.

### Methodology

- Original model reconstructed based on the published paper and public datasets when available.
- PLDR implemented by partitioning the uncertainty space and associating linear responses within each partition.
- Comparative experiments performed on modified IEEE 300-bus scenarios using Julia + JuMP + Gurobi.

### Results Summary

- The PLDR approach showed promising improvements in adaptability, especially in highly volatile or nonlinear uncertainty regimes.
- However, it introduced significant computational challenges, particularly in larger systems or with fine-grained partitions.

## Academic Integrity Notice

This repository is made available for educational and academic dissemination only.  
**Any direct use of this material for coursework or derivative submissions without proper citation is strictly forbidden.**  
All work presented here was developed independently and made public after course completion.

## Acknowledgment

This work was inspired by the methodology of Rodrigues, Street, and Arroyo (2024), and builds upon their contributions to multistage stochastic programming in power systems.

