SloppyCell Example
==================

In this example, we'll first reproduce the results of and build upon the work in avaiable for free `on Pubmed`__ 

.. _Tyson1991: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=1831270&query_hl=1
__ Tyson1991

Creating the Networks
---------------------

Import all the SloppyCell tools::

  from SloppyCell.ReactionNetworks import *


Load up the basic model from an SBML file obtained from `BioModels.net`_. ``base_net`` is the Python variable name for our network, and ``base`` is the id we give our Network for later incorporation into a Model.::

  base_net = IO.from_SBML_file('BIOMD0000000005.xml', 'base')

The SBML file is a little sloppy. Let's clean it up a bit.

  The variables ``cell`` and ``EmptySet`` don't change, but aren't explicitly
  declared constant in the SBML file. We'll set them constant here to slightly
  speed up our integrations.::

    base_net.set_var_constant('cell', True)
    base_net.set_var_constant('EmptySet', True)

  The total amount of cyclin, denoted ``YT``, appears in several of the paper's
  plots. Let's define such a species::

    base_net.add_species('YT', 'cell', name = 'total cyclin')
    base_net.addAssignmentRule('YT', 'Y+YP+pM+M')

  Later we'll be running some optimization routines on our Model. Two parameters
  in the original Tyson network are 0. Let's set those parmeters to not be 
  optimizable.

    base_net.set_var_optimizable('k5notP', False)
    base_net.set_var_optimizable('k2', False)

Now that we've made our changes, we need to ``compile`` the network. This will
parse the components of the model and generate the various differential 
equations we'll use::
  
  base_net.compile()

Let's output our equations to LaTeX form. The results of this can be found `here <TeX/eqns.pdf>`_::

  IO.eqns_TeX_file(base_net, 'TeX/eqns.tex')

.. figure:: 3.png
   :scale: 50
   :align: center

   Our reproduction of figure 3 from Tyson1991_.

.. _BioModels.net: http://www.biomodels.net
