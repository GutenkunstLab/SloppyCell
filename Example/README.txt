SloppyCell Example
==================

In this example, we'll first reproduce the results of and build upon the work in avaiable for free `on Pubmed`__ 

.. _Tyson1991: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=1831270&query_hl=1
__ Tyson1991

Creating the Networks
---------------------

Import all the SloppyCell tools::

  from SloppyCell.ReactionNetworks import *


Load up the basic model from an SBML file obtained from `BioModels.net`_. ``base`` is the Python variable name for our network, and ``base`` is the id we give our Network for later incorporation into a Model.::

  base_net = IO.from_SBML_file('BIOMD0000000005.xml', 'base')

The SBML file is a little sloppy. Let's clean it up a bit.

  The variables ``cell`` and ``EmptySet`` don't change, but aren't explicitly
  declared constant in the SBML file. We'll set them constant here to slightly
  speed up our integrations.::

    base_net.set_var_constant('cell', True)
    base_net.set_var_constant('EmptySet', True)

  The total amount of cyclin, denoted ``YT``, appears in several of the paper's
  plots. Let's define such a species. First, we create the species ``YT``, which
  exists in the ``cell`` compartment. We also give it a more descriptive name 
  for plots and TeX output.::

    base_net.add_species('YT', 'cell', name = 'total cyclin')
    base_net.add_assignment_rule('YT', 'Y+YP+pM+M')

  Later we'll be running some optimization routines on our Model. Two parameters
  in the original Tyson network are 0. Let's set those parameters to not be 
  optimizable.

    base_net.set_var_optimizable('k5notP', False)
    base_net.set_var_optimizable('k2', False)

Let's output our equations to LaTeX form. The results can be found 
`here <TeX/eqns.pdf>`_::

  IO.eqns_TeX_file(base_net, 'TeX/eqns.tex')

Now that we have a well-constructed base network, we can create some of the 
modified networks used in the paper. In figure 3b, there's an example of a
network being perturbed. 

  We copy our original network, giving it a new name in the process.::

    perturbed_net = base_net.copy('perturbed')

  We want to set ``k6`` to 2.0. Let's also make sure that it doesn't change
  during optimization.::

    perturbed_net.set_initial_var_value('k6', 2.0)
    perturbed_net.set_var_optimizable('k6', False)

  Figure 3b starts from a fixed point, so we'll find that. We'll integrate
  our equations for a while, then run a Newton algorithm to find the fixed
  point.::

    Dynamics.integrate(perturbed_net, [0, 100])
    Dynamics.dyn_var_fixed_point(perturbed_net)

  The we set the initial conditions of our dynamic variables to be equal to the
  fixed point.::

    perturbed_net.set_dyn_var_ics(fp)

  Finally, we add a number of events to perturb our network. In each case, the
  first argument is an identifier. The second is a trigger, when this expression
  becomes True, the event will fire. For example, the first event fires when
  the time becomes greater than 10. The final argument is a dictionary of
  event assignments, giving the new values the variables get when the event
  executes. In this case, ``M`` is increased by various factors.::

    perturbed_net.add_event('perturb1', 'gt(time, 10)', {'M': '1.2*M'})
    perturbed_net.add_event('perturb2', 'gt(time, 30)', {'M': '2*M'})
    perturbed_net.add_event('perturb3', 'gt(time, 60)', {'M': '2.4*M'})

  We add a parameter to the network, with the value 116. By default, parameters
  are constant and optimizable.::

    growth_net.add_parameter('Td', 116)

  We add a rate rule for ``k6``. This sets the derivative wrt to time of the
  first argument to be equal to the second. Here we have exponential decay.::

    growth_net.add_rate_rule('k6', '-0.693 * k6/Td')
    growth_net.set_initial_var_value('k6', 2.5)

  And we add an event. The trigger here is when ``k6`` becomes less than 1.3.::

    growth_net.add_event('gene_replication', 'lt(k6, 1.3)', {'k6': '2 * k6'})

.. figure:: 3.png
   :scale: 50
   :align: center

   Our reproduction of figure 3 from Tyson1991_.

.. _BioModels.net: http://www.biomodels.net
