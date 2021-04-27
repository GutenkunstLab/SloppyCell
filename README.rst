SloppyCell
==========

SloppyCell is a software environment for simulation and analysis of biomolecular networks. A particular strength of SloppyCell is estimating parameters by fitting experimental data and then calculating the resulting uncertainties on parameter values and model predictions.

SloppyCell was initially developed by Ryan Gutenkunst in the research group of `Jim Sethna <http://sethna.lassp.cornell.edu/>`_. Ryan now has his `own group <http://gutengroup.mcb.arizona.edu>`_, where SloppyCell continues to be developed. Most recently, an XML interface was been developed in collaboration with `Michael Blinov <http://apache.cam.uchc.edu/mblinov/>`_, to facilitate integration with other modeling frameworks, such as `The Virtual Cell <http://vcell.org>`_.

Note
----
Recently, there appears to have been a change to Python distutils such that networks cannot be compiled in an imported module. If attempt is made to compile a network in an imported module, SloppyCell will hang.

Features
--------
* support for much of the Systems Biology Markup Language (SBML) level 2 version 3
* deterministic and stochastic dynamical simulations
* sensitivity analysis without finite-difference derviatives
* optimization methods to fit parameters to experimental data
* simulation of multiple related networks sharing common parameters
* stochastic Bayesian analysis of parameter space to estimate uncertainties associated with optimal fits
