#
# This script makes the htdocs directory for SloppyCell.
#
# It downloads the cvs source, builds it, then runs the tests.
# Then it builds the documentation.
#

import os
import glob
import shutil
import time

# Name of the temporary directory to export SloppyCell to
TEMP='.tmp_SloppyCell'

# Change to the HOME directory
home = os.getenv('HOME')
os.chdir(home)

# Delete the temporary directory, remake it, then change into it.
os.system('rm -rf %s' % TEMP)
os.mkdir(TEMP)
os.chdir(TEMP)

# Export the code from CVS
os.system('cvs -d:pserver:anonymous@sloppycell.cvs.sourceforge.net:/cvsroot/sloppycell login')
os.system('cvs -z3 -d:pserver:anonymous@sloppycell.cvs.sourceforge.net:/cvsroot/sloppycell co -P SloppyCell')

# Make this directory the first in the PYTHONPATH
old_python_path = os.getenv('PYTHONPATH')
new_python_path = '%s:%s' % (os.path.abspath(os.curdir), old_python_path)
os.putenv('PYTHONPATH', new_python_path)

# Enter the SloppyCell directory
os.chdir('SloppyCell')

# Build SloppyCell and install the libraries to their proper places
build_failure = os.system('python setup.py build >& build.output')
if build_failure:
    os.sys.exit()
os.system('python setup.py install --install-lib=..')

# Run test.py
os.chdir('Testing')
test_failure = os.system('python test.py --debugSC=test.ouput')
if test_failure:
    os.sys.exit()

os.chdir('..')

# Done running SloppyCell, restore PYTHONPATH
os.putenv('PYTHONPATH', old_python_path)

#
# On to documentation...
#

# Change the version info
os.system(r"""cat __init__.py | sed "s/CVS/%s/g" > __temp_dork.py"""
          % time.asctime().strip())
shutil.move('__temp_dork.py', '__init__.py')

os.system(r"""cat Doc/index.html | sed "s/CVS_GEN_DATE/%s/g" > Doc/.tmp1"""
          % time.asctime())
os.system(r"""cat Doc/.tmp1 | sed "s/API_DOC_GEN_DATE/%s/g" > Doc/.tmp2"""
          % time.asctime())
shutil.move('Doc/.tmp2', 'Doc/index.html')

# Make the source distribtions
os.system('python setup.py sdist --formats=gztar,zip')

os.mkdir('htdocs')
shutil.move('dist/SloppyCell-CVS.tar.gz', 'htdocs/SloppyCell-CVS.tar.gz')
shutil.move('dist/SloppyCell-CVS.zip', 'htdocs/SloppyCell-CVS.zip')

os.chdir('..')
os.system('epydoc SloppyCell/ -v -o SloppyCell/htdocs/api --html --no-frames --url=http://sloppycell.sf.net --name=SloppyCell --docformat=plaintext --parse-only --exclude=Testing* --exclude=Example* --show-imports --include-log')

os.chdir('SloppyCell/Doc/devel')
os.system('latex devel')
os.system('bibtex devel')
os.system('latex devel')
os.system('latex devel')
os.system('dvipdf devel')

os.chdir('../user')
os.system('latex user')
os.system('bibtex user')
os.system('latex user')
os.system('latex user')
os.system('dvipdf user')

os.chdir('../..')
os.system('cp Doc/index.html Doc/user/user.pdf Doc/devel/devel.pdf htdocs/.')

os.system('chmod -R a+rx htdocs')

# Use this command to copy things up to SourceForge.
upload_command = 'scp -r $HOME/%s/SloppyCell/htdocs jepetto@shell.sf.net:/home/groups/s/sl/sloppycell/.' % TEMP
print 'To upload the site to SourceForge, use scp like: %s' % upload_command
print 'Use your username in place of jepetto.'
