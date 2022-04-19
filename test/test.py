from __future__ import print_function
import glob
import unittest
import sys
from SloppyCell import ReactionNetworks
from SloppyCell import Utility
from SloppyCell import disable_c

from SloppyCell.ReactionNetworks  import Network_mod 

def run_all_tests():
    all_tests = unittest.TestSuite()

    testfiles = glob.glob('test_*.py')
    mesgs = []
    all_test_mods = []
    for file in testfiles:
        module = file[:-3]
        mod = __import__(module)
        all_test_mods.append(mod)
        print("all tests", all_test_mods)
        if hasattr(mod, 'suite'):
            all_tests.addTest(mod.suite)
    print(all_tests)

    if not '-v' in sys.argv:
        Utility.disable_warnings()
    if not disable_c:
        print('*' * 80)
        print('Running tests with C compilation enabled.')
        print('*' * 80)
        unittest.TextTestRunner(verbosity=2).run(all_tests)
    ReactionNetworks.Network_mod.Network.disable_c = True
    print('*' * 80)
    print('Running tests with C compilation disabled.')
    print('*' * 80)
    unittest.TextTestRunner(verbosity=2).run(all_tests)

    for mod in all_test_mods:
        if hasattr(mod, 'message'):
            print()
            print(mod.message)

if __name__ == '__main__':
    run_all_tests()