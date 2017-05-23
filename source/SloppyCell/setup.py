# This gets f2py_signatures installed in the correct place. See 
# http://groups.google.com/group/comp.lang.python/browse_thread/thread/35ec7b2fed36eaec/2105ee4d9e8042cb
from distutils.command.install import INSTALL_SCHEMES
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib'] 

# Importing these adds a 'bdist_mpkg' option that allows building binary packages on OS X.
try:
    import setuptools
    import bdist_mpkg
except ImportError:
    pass


import scipy
if hasattr(scipy, 'Numeric'):
    # Using old scipy
    import scipy_distutils.core as core
    data_files = [('SloppyCell/ReactionNetworks', 
                   ['ReactionNetworks/f2py_signatures.pyf',
                    'ReactionNetworks/f2py_signatures_no_derivs.pyf',
                    'ReactionNetworks/mtrand.h',
                    'ReactionNetworks/mtrand.c'])]
else:
    # Using new scipy
    import numpy.distutils.core as core

    from numpy.distutils.misc_util import Configuration
    # Annoyingly, distutils doesn't have an option to recursively include
    # a data directory. We can use numpy's Configuration class, however,
    # to do the recursion for us.
    config = Configuration('SloppyCell')
    config.add_data_files('ReactionNetworks/f2py_signatures.pyf',
                          'ReactionNetworks/f2py_signatures_no_derivs.pyf',
                          'ReactionNetworks/mtrand.h',
                          'ReactionNetworks/mtrand.c')
    config.add_data_dir('Doc')
    config.add_data_dir('ddaskr')
    config.add_data_dir('Example')
    data_files = config.todict()['data_files']

# This is a kludge so that sdist works on AFS systems.
# AFS doesn't allow hard-linking across directories, but if we're on linux
#  we still have a os.link function (which is all distutils check for).
# If we delete os.link, sdist will copy rather than link.
# See http://mail.python.org/pipermail/distutils-sig/2002-October/002990.html
#  for more detail, and a patch that apparently never got applied to trunk
#  distutils
# I don't know why this kludge wasn't needed for Python 2.3
import os
if hasattr(os, 'link'):
    del os.link

# These packages include some BLAS functions in them. Strangely, at least
#  on CCMR, it seems (very slightly), faster to use the included ones rather
#  than linking against LAPACK.
daskr = core.Extension(name = 'SloppyCell._daskr',
                       sources = ['daskr.pyf', 'ddaskr/ddaskr.f', 
                                  'ddaskr/daux.f', 'ddaskr/dlinpk.f'])

misc_c = core.Extension(name = 'SloppyCell.misc_c',
                        sources = ['misc_c.c', 'misc_c.pyf'])

core.setup(name='SloppyCell',
           version='CVS',
           author='Ryan Gutenkunst',
           author_email='rng7@cornell.edu',
           url='http://sloppycell.sourceforge.net',
           packages=['SloppyCell', 
                     'SloppyCell.ReactionNetworks', 
                     'SloppyCell.Testing',
                     'SloppyCell.ExprManip',
                     'SloppyCell.Vandermonde',
                     ],
           package_dir={'SloppyCell': ''},
           data_files=data_files,

           ext_modules = [daskr, misc_c]
           )
