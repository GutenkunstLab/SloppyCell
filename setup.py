import glob
from scipy_distutils.core import setup, Extension

# We want the data files for the various xsl transforms included in the
# installation tree. That keeps the directory structure clean, and they're
# easy to find then.
#
# The ideal way to do this would be to use the package_data argument to setup,
# but that's not implemented in 2.3. For now, we'll use this hack from Pete
# Shinners.
from distutils.command.install_data import install_data
class smart_install_data(install_data):
    def run(self):
        #need to change self.install_dir to the library dir
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)

lsodar = Extension(name = 'SloppyCell._lsodar',
                   sources = ['lsodar.pyf', 'odepack/opkdmain.f', 
                              'odepack/opkda1.f', 'odepack/opkda2.f'])

setup(name='SloppyCell',
      version='0.2',
      author='Ryan Gutenkunst',
      author_email='rng7@cornell.edu',
      url='http://sloppycell.sourceforge.net',
      packages=['SloppyCell', 
                'SloppyCell.ReactionNetworks', 
                'SloppyCell.Testing'
                ],
      package_dir={'SloppyCell': ''},

      # More of the package_data hack
      cmdclass = {'install_data': smart_install_data},
      data_files=[('SloppyCell/ReactionNetworks', 
                   ['ReactionNetworks/sbml_l2v1_todot.xsl']),
                  ('SloppyCell/ReactionNetworks/xsltml_2.0', 
                   glob.glob('ReactionNetworks/xsltml_2.0/*.xsl') +\
                   ['ReactionNetworks/xsltml_2.0/README']),
                  ('SloppyCell/Doc', glob.glob('Doc/*'))],

      # The 2.4 way to replace the hack
      #package_data={'SloppyCell.ReactionNetworks': ['sbml_l2v1_todot.xsl', 
      #                                              'xsltml_2.0/*']},
      ext_modules = [lsodar]
      )
