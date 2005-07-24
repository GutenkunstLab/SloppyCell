import glob
import unittest

def run_all_tests():
    all_tests = unittest.TestSuite()

    testfiles = glob.glob('test_*.py')
    mesgs = []
    for file in testfiles:
        module = file[:-3]
        mod = __import__(file[:-3])
        if hasattr(mod, 'suite'):
            all_tests.addTest(mod.suite)
        if hasattr(mod, 'message'):
            mesgs.append(mod.message)

    unittest.TextTestRunner(verbosity=2).run(all_tests)

    print
    for mesg in mesgs:
        print mesg
        print

if __name__ == '__main__':
    run_all_tests()
