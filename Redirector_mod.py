import os
import sys
import tempfile

# Adapted from code by Robert Kern. 
# See http://osdir.com/ml/python.ipython.devel/2005-08/msg00014.html

STDOUT = 1
STDERR = 2

class Redirector(object):
    def __init__(self, fd=STDOUT):
        self.fd = fd
        self.started = False

    def start(self):
        if not self.started:
            self.flush()
            self.tmpfd, self.tmpfn = tempfile.mkstemp()

            self.oldhandle = os.dup(self.fd)
            os.dup2(self.tmpfd, self.fd)
            os.close(self.tmpfd)

            self.started = True

    def flush(self):
        if self.fd == STDOUT:
            sys.stdout.flush()
        elif self.fd == STDERR:
            sys.stderr.flush()

    def stop(self):
        if self.started:
            self.flush()
            os.dup2(self.oldhandle, self.fd)
            os.close(self.oldhandle)
            tmpr = open(self.tmpfn, 'rb')
            output = tmpr.read()
            tmpr.close()  # this also closes self.tmpfd
            try:
                os.unlink(self.tmpfn)
            except OSError:
                pass

            self.started = False
            return output
        else:
            return None

# Adapted from http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
# Added since the above redirector doesn't seem to
# work in ipython notebook for Network.fun_f2py.
class hideStdout(object):
    def __init__(self):
        self.devnull = open(os.devnull, 'w')
        self._stdout = self.devnull
    
    def start(self):
        self.old_stdout = sys.stdout
        self.old_stdout.flush()
        sys.stdout = self._stdout
    
    def stop(self):
        self._stdout.flush()
        sys.stdout = self.old_stdout
