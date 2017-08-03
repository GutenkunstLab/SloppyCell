# Importing these adds a 'bdist_mpkg' option that allows building binary packages on OS X.
try:
    import setuptools
    import bdist_mpkg
except ImportError:
    pass

import numpy.distutils.core as core

daskr = core.Extension(name = 'SloppyCell._daskr',
                       sources = ['SloppyCell/daskr.pyf',
                                  'ddaskr/ddaskr.c',
                                  'ddaskr/daux.c',
                                  'ddaskr/dlinpk.c',
                                  'ddaskr/ddaskr_types.h'])

misc_c = core.Extension(name = 'SloppyCell.misc_c',
                        sources = ['SloppyCell/misc_c.c',
                                   'SloppyCell/misc_c.pyf'])

rxn_data_files = ['SloppyCell/ReactionNetworks/{0}'.format(_) for _
                  in ['mtrand.c', 'mtrand.h', 'f2py_signatures.pyf',
                      'f2py_signatures_no_derivs.pyf']]

core.setup(name='SloppyCell',
           version='1.1.0',
           author='Ryan Gutenkunst',
           author_email='rgutenk@email.arizona.edu',
           url='http://sloppycell.sourceforge.net',
           packages=['SloppyCell', 'SloppyCell.ExprManip',
                     'SloppyCell.ReactionNetworks',
                     'SloppyCell.Vandermonde'],
           ext_modules = [daskr, misc_c],
           zip_safe=False,
           data_files=[('SloppyCell/ReactionNetworks', rxn_data_files)],
           include_package_data = True,
          )
