import setuptools
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

# Collect files needed for dynamic compilation
rxn_data_files = ['SloppyCell/ReactionNetworks/{0}'.format(_) for _
                  in ['mtrand.c', 'mtrand.h', 'f2py_signatures.pyf',
                      'f2py_signatures_no_derivs.pyf']]

# For reading a file into lon_description 
import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

core.setup(name='SloppyCell',
           version='1.1.0.dev1',
           author='Ryan Gutenkunst',
           author_email='rgutenk@email.arizona.edu',
           url='http://sloppycell.sourceforge.net',
           license='BSD',
           long_description=read('README.rst'),
           classifiers=['Development Status :: 5 - Production/Stable',
                        'Environment :: Console',
                        'Intended Audience :: Science/Research',
                        'License :: OSI Approved :: BSD License',
                        'Operating System :: OS Independent',
                        'Programming Language :: Python :: 3',
                        'Programming Language :: C',
                        'Topic :: Scientific/Engineering :: Bio-Informatics'
                       ],
           scripts=['SloppyCell/sloppycell.py'],
           keywords='systems biology, sloppy modeling',
           requires=['scipy', 'numpy', 'matplotlib'],
           install_requires=['numpy', 'scipy'],
           packages=['SloppyCell', 'SloppyCell.ExprManip',
                     'SloppyCell.ReactionNetworks',
                     'SloppyCell.Vandermonde'],
           ext_modules = [daskr, misc_c],
           # Ensure files needed for dynamic compilation are installed
           zip_safe=False,
           data_files=[('SloppyCell/ReactionNetworks', rxn_data_files)],
           include_package_data = True,
          )
