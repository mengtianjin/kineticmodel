\# CO2RR Kinetic Modeling



This repository contains a MATLAB-based kinetic modeling workflow for simulating and optimizing the electrocatalytic reduction of CO2 (CO2RR). It fits experimental partial current densities (CO and H2) using a microkinetic model, extracting intrinsic rate constants and transfer coefficients.



\## File Structure



\* `main\_kinetic\_modeling.m`: The primary script. It handles data initialization, configures the `MultiStart` global optimization environment, runs the solver, and calculates the final simulated coverages and current densities.

\* `ObjectiveFun\_Log.m`: The custom objective function. It calculates the steady-state surface coverages ($\\theta\_{CO}$, $\\theta\_{H}$, $\\theta\_{bare}$) and computes the log-space residuals between the experimental and modeled reaction rates.



\## Requirements



\* MATLAB (R2020a or newer recommended)

\* Optimization Toolbox (for `lsqnonlin`)

\* Global Optimization Toolbox (for `MultiStart`)

\* Parallel Computing Toolbox (optional but highly recommended for speeding up the multi-start process)





