SloppyCell is a software environment for simulation and analysis of 
biomolecular networks.  SloppyCell is written in Python, and leverages 
efficient numerical algorithms provided by the SciPy package (www.scipy.org).
Features of SloppyCell include:

- deterministic and stochastic dynamical simulations
- simulation of multiple related networks sharing common values of 
  parameters
- support for the Systems Biology Markup Language (SBML) level 2 version 1
- forward and backward adjoint sensitivity analysis
- parameter optimization methods to fit parameters to experimental data
- stochastic Bayesian analysis of parameter space to estimate error bars
  associated with optimal fits

In addition to SciPy, SloppyCell uses pylab (matplotlib.sourceforge.net) for
graphics, and libsbml (www.sbml.org) for SBML interfaces.

The top level of the SloppyCell package provides support for describing
Models, Calculations and Experiments.  Models combine Calculations and 
Experiments to support the definition of a cost function (the residual
misfit between Calculations and Experiments).  The most widely used
type of Calculation is that described by a ReactionNetwork.  The 
ReactionNetworks subdirectory contains code to support the construction
of networks from sets of chemicals, reactions, and parameters.



