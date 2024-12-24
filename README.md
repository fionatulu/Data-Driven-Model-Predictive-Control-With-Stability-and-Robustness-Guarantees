# Data-Driven-Model-Predictive-Control-With-Stability-and-Robustness-Guarantees
This repository contains the implementation and reproduction of the findings from the paper titled "Data-Driven Model Predictive Control With Stability and Robustness Guarantees". The goal of this project is to verify the results and methodology presented in the original research.

According to the paper, with the trajectory of collected input-output data, along with the upper bound on the system's order, future input-output data can be calculated with Hankel matrix. This work also proved the recursive feasibility and exponential stability of the method.

## Code Reference
The code in this repository is heavily inspired by the work of [shiivashaakeri], which can be found at [(https://github.com/shiivashaakeri/Data-Driven-Model-Predictive-Control-MPC-with-Stability-and-Robustness-Guarantees)]. We have made modifications according to the original paper of this project and to ensure compatibility with Matlab.

## Implementation Files
This repository contains the following matlab scripts:
- henkel.m: Module for handling Hankel matrix operations.
- MDL_sim_prestab.m: Script to simulate the model with pre-stabilization and get the priori measured data trajectory of length N. 
- robust_mpc.m: The main script implementing the robust data-driven MPC algorithm.
Moreover, the experiment is applied to the four tank system, the system file is set as defined in the original paper. The definition file is:
- tank_sys2.mat: four tank system file.
