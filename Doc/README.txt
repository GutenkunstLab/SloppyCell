SloppyCell
==========

SloppyCell is a software environment for simulation and analysis of biomolecular
networks developed by the groups of `Jim Sethna`_ and `Chris Myers`_ at `Cornell University`_.

Examples of models developed in SloppyCell can be found at `Jim Sethna's Gene Dynamics page <http://www.lassp.cornell.edu/sethna/GeneDynamics/>`_.

.. topic:: News

   January 23rd, 2007: **SloppyCell 0.6 released!**

      Lots of work over the past few months. The parallel processing code has
      been completely rewritten. Preliminary support for algebraic rules is also
      in (many thanks to Jordan Atlas) for his help.

      Grab it from our `download page <http://sourceforge.net/project/showfiles.php?group_id=140498>`_.

   September 22nd, 2006: **New windows installer**

      An installer for Python 2.4 on Windows is now available. See the 
      `installation guide <INSTALL.html>`_ for details.

   September 15th, 2006: **Speed-up for new numpy**

      The latest development versions of numpy have (at my request) added a
      method to speed up array access. In our case, this can speed up 
      integrations by ~30%. The latest CVS version of SloppyCell supports this
      new interface while also supporting the older SciPy/numpy.

   September 13th, 2006: **"New" SciPy compatibility**

      The CVS version of SloppyCell has been updated to work with new
      versions of numpy and SciPy. Testing has been on (SciPy 0.5.1, 
      numpy 1.0b5, matplotlib 0.87.5) and on (SciPy 0.3.2, Numeric 23.1, 
      matplotlib 0.83.2).


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
- `Usage example <Example/README.html>`_ (incomplete, see the Example directory 
  in the source)
- `Mailing list`_
- Sourceforge `project page`_
- CVS snapshot (`.tar.gz <SloppyCell-CVS.tar.gz>`_, `.zip <SloppyCell-CVS.zip>`_) generated CVS_GEN_DATE
- `API documentation <api/index.html>`_ generated API_DOC_GEN_DATE

.. topic:: Citing SloppyCell:

   There are no papers yet on SloppyCell itself. For now, please reference this website:

      Ryan N. Gutenkunst, Fergal P. Casey, Joshua J. Waterfall, Jordan C. Atlas, Robert S. Kuczenski, Chris R. Myers, and James P. Sethna. SloppyCell. http://sloppycell.sourceforge.net/ (2007)
  

$Date$

.. _Cornell University: http://www.cornell.edu
.. _Mailing list: http://lists.sourceforge.net/lists/listinfo/sloppycell-users
.. _changelog: changelog.html
.. _project page: http://www.sourceforge.net/projects/sloppycell/
.. _Jim Sethna: http://www.lassp.cornell.edu/sethna/
.. _Chris Myers: http://www.tc.cornell.edu/~myers/
.. _sloppycell-users: mailto:sloppycell-users@lists.sourceforge.net
.. _Ryan Gutenkunst: mailto:rng7@cornell.edu
