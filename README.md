# integrate-hse
Integrates the hydrostatic equilibrium equations of stellar structure. A polytropic EOS relates density and pressure.
Currently configured for neutron star interiors.

Uses the SUNDIALS CVODE package to apply the Adams-Moulton linear multistep method for the non-stiff ODE system.
At each step, the nonlinear solution equations are solved via Newton iteration using the analytic Jacobian matrix.

Originally written for Mike Zingale's Stars course at Stony Brook University: 
http://bender.astro.sunysb.edu/classes/stars/
