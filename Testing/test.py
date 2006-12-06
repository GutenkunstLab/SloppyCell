import glob
import unittest

def run_all_tests():
    all_tests = unittest.TestSuite()

    testfiles = glob.glob('test_*.py')
    mesgs = []
    all_test_mods = []
    for file in testfiles:
        module = file[:-3]
        mod = __import__(module)
        all_test_mods.append(mod)
        if hasattr(mod, 'suite'):
            all_tests.addTest(mod.suite)

    unittest.TextTestRunner(verbosity=2).run(all_tests)

    for mod in all_test_mods:
        if hasattr(mod, 'message'):
            print
            print mod.message

if __name__ == '__main__':
    run_all_tests()
