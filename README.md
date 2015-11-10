# integrate-hse
Integrates the hydrostatic equilibrium equations of stellar structure. A polytropic EOS relates density and pressure.
Currently configured for neutron star interiors.

Uses the SUNDIALS CVODE package to apply the Adams-Moulton linear multistep method for the non-stiff ODE system.
At each step, the nonlinear solution equations are solved via Newton iteration using the analytic Jacobian matrix.

Originally written for Mike Zingale's Stars course at Stony Brook University: 
http://bender.astro.sunysb.edu/classes/stars/

To generate the Mass-Radius and Mass-Central Density curves:
1) make
2) python run_m-r.py
3) python plot_m-r.py

Compiles with gfortran 5.1.1 and uses SUNDIALS 2.6.2, LAPACK 3.5.0, & BLAS 3.5.0.
Plotting scripts use Python 2.7.10 provided by Anaconda 2.3.0.