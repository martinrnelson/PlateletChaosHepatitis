This repository contains codes related to the paper "Platelet-driven routes to chaos in a model of hepatitis", by M.R. Nelson, J.M. Gibbins and J.L. Dunster, submitted to Chaos, Solitons and Fractals.

The key files to run in order to reproduce the Figures from this manuscript are explained below. The remaining Matlab files in the repository are called by the files below, and brief documentation can be found within these files.

XPP FILES:

- PlateletChaosHepatitis_oneRhoParameter.ode: implementation of the ODE system (1) in the manuscript.
- PlateletChaosHepatitis_manyRhoParameters.ode: implementation of the ODE system (3) in the manuscript.

MATLAB FILES:

For the "single rho" model of Section 4:

- odeSolver_singleRho.m: Numerical integrator (produces e.g. Figure 3).
- plotMaxMinScatters.m: Produces Figure 4(c,d).
- plotAttractor.m: Produces e.g. Figures 5 and 6.
- calculateLLE_singleRho.m: Calculates and plots Largest Lyapunov Exponents; produces e.g. Figure 7(d,e,f).

For the "multiple rhos" model of Section 5:

- mcmc_manyRhos.m: Run a Metropolis-Hastings algorithm to explore 7D rho-space seeking maximal LLE; produces data for the upper-right panels of Figure 9. (See output file MCMC_results.mat.)
- calculateLLE_manyRhos_maxLLE.m: Take the outputs from mcmc_manyRhos.m, take the point of maximal LLE of those sampled and compute the local LLE landscape on 2D slices of rho-space; produces data for the lower-left panels of Figure 9. (See output file lyapunovCalc_maxMCMC.mat.)
- plotManyRhosFigures.m: Plot the results of both of the above codes to produce Figure 9.

Any queries regarding the use of these functions can be sent to martin.nelson@ntu.ac.uk at any time.


