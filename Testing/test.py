import glob
import unittest

import logging
import sys
if '--debug' in sys.argv:
    logging.getLogger().setLevel(logging.DEBUG)

def run_all_tests():
    all_tests = unittest.TestSuite()

    testfiles = glob.glob('test_*.py')
    mesgs = []
    all_test_mods = []
    for file in testfiles:
        module = file[:-3]
        mod = __import__(file[:-3])
        if hasattr(mod, 'suite'):
            all_test_mods.append(mod)
            all_tests.addTest(mod.suite)

    unittest.TextTestRunner(verbosity=2).run(all_tests)

    for mod in all_test_mods:
        if hasattr(mod, 'message'):
            mesgs.append(mod.message)

    print
    for mesg in mesgs:
        print mesg
        print

if __name__ == '__main__':
    run_all_tests()
