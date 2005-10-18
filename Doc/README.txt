SloppyCell
==========

SloppyCell is a software environment for simulation and analysis of biomolecular
networks developed by the groups of `Jim Sethna`_ and `Chris Myers`_ at `Cornell University`_.

Examples of models developed in SloppyCell can be found at `Jim Sethna's Gene Dynamics page <http://www.lassp.cornell.edu/sethna/GeneDynamics/>`_.

.. topic:: News

   October 18th, 2005: **Welcome ICSB Participants!**

      Perhaps we have some visitors who saw our posters at the 6th International
      Conference on Systems Biology. If you're interested in our methods, please
      download the CVS tarball below and check it out.

      Unfortunately, we weren't able to get all we wanted finished in
      time for the conference. More complete documentation is forthcoming, as
      is an updated Windows build. In the meantime, you can look in
      SloppyCell/Example for a script that uses almost all of our functionality.
      
      Please email the sloppycell-users `mailing list`_ for assistance, 
      or contact `Ryan Gutenkunst`_ directly.

   July 4th, 2005: **SloppyCell 0.1 Released**

      This is the initial public release of SloppyCell. Please see the 
      changelog_ for changes from recent development versions.


.. topic:: Features:

   - support for much of the `Systems Biology Markup Language (SBML) 
     <http://sbml.org>`_ level 2 version 1 (`details <SBML_support.html>`_)
   - deterministic dynamical simulations (stochastic on the way)
   - sensitivity analysis
   - optimization methods to fit parameters to experimental data
   - simulation of multiple related networks sharing common parameters
   - stochastic Bayesian analysis of parameter space to estimate error bars
     associated with optimal fits

.. The top level of the SloppyCell package provides support for describing Models, Calculations and Experiments.  Models combine Calculations and Experiments to support the definition of a cost function (the residual misfit between Calculations and Experiments).  The most widely used type of Calculation is that described by a ReactionNetwork.  The ReactionNetworks subdirectory contains code to support the construction of networks from sets of chemicals, reactions, and parameters.

- `Installation guide <INSTALL.html>`_
- `Usage example <Example/README.html>`_
- `Mailing list`_
- Sourceforge `project page`_
- `CVS snapshot <SloppyCell-CVS.tar.gz>`_ generated nightly
- `API documentation <api/index.html>`_ generated nightly

.. topic:: Citing SloppyCell:

    There are no papers yet on SloppyCell itself. For now, please reference this website:

    Ryan Gutenkunst, Fergal Casey, Josh Waterfall, Kevin Brown, Chris Myers, and Jim Sethna. SloppyCell. http://sloppycell.soureforge.net/ (2005)
  

$Date$

.. _Cornell University: http://www.cornell.edu
.. _Mailing list: http://lists.sourceforge.net/lists/listinfo/sloppycell-users
.. _changelog: changelog.html
.. _project page: http://www.sourceforge.net/projects/sloppycell/
.. _Jim Sethna: http://www.lassp.cornell.edu/sethna/
.. _Chris Myers: http://www.tc.cornell.edu/~myers/
.. _sloppycell-users: mailto:sloppycell-users@lists.sourceforge.net
.. _Ryan Gutenkunst: mailto:rng7@cornell.edu
