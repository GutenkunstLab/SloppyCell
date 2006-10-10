try:
    import numpy.distutils.core as core
except ImportError:
    import scipy_distutils.core as core

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

lsodar = core.Extension(name = 'SloppyCell._lsodar',
                        sources = ['lsodar.pyf', 'odepack/opkdmain.f', 
                                   'odepack/opkda1.f', 'odepack/opkda2.f'])

core.setup(name='SloppyCell',
           version='CVS',
           author='Ryan Gutenkunst',
           author_email='rng7@cornell.edu',
           url='http://sloppycell.sourceforge.net',
           packages=['SloppyCell', 
                     'SloppyCell.ReactionNetworks', 
                     'SloppyCell.Testing',
                     'SloppyCell.ExprManip'
                     ],
           package_dir={'SloppyCell': ''},

           ext_modules = [lsodar]
           )
