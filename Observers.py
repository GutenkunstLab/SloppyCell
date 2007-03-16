import os
import time

import scipy

import Utility

def print_hess_elements(**args):
    if args['event'] == 'hessian element':
        elem = args['element']
        ii = args['i']
        jj = args['j']

        print 'hessian element %i, %i: %g' % (ii, jj, elem)

class CostPrinter:
    def __init__(self, skip=1, print_params=False, print_best_params=False):
        self.skip = skip
        self.print_params = print_params
        self.print_best_params = print_best_params
        self.reset()

    def __call__(self, **args):
        if args['event'] == 'evaluation':
            cost = args['cost']
            params = args['params']

            if cost < self.lowest_cost:
                self.lowest_cost = cost
                self.best_params = params.copy()
            if self.ii % self.skip == 0:
                print 'call %i: cost: %g, best so far: %g' % (self.ii, cost, 
                                                              self.lowest_cost)
                os.sys.stdout.flush()
                if self.print_params:
                    print params
                if self.print_best_params:
                    print self.best_params
            self.ii += 1

    def reset(self):
        self.ii = 0
        self.lowest_cost = scipy.inf
        self.best_params = None

def print_all_costs(**args):
    if args['event'] == 'cost':
        print args['cost']

class CostEmailer:
    def __init__(self, interval, from_addr, to_addr):
        self.interval = interval
        self.from_addr, self.to_addr = from_addr, to_addr
        self.reset()

    def __call__(self, **args):
        if args['event'] == 'evaluation':
            cost = args['cost']
            params = args['params']

            if cost < self.lowest_cost:
                self.lowest_cost = cost
                self.best_params = params.copy()

            if time.time() - self.last_sent > self.interval * 3600:
                lines = []
                lines.append('Best cost so far: %f' % self.lowest_cost)
                lines.append('Corresponding to parameters: %s' 
                             % str(self.best_params))
                msg = os.linesep.join(lines)

                Utility.send_email(self.to_addr, self.from_addr,
                                   "SloppyCell job update",
                                   msg)
                self.last_sent = time.time()

            self.ii += 1

    def reset(self):
        self.ii = 0
        self.lowest_cost = scipy.inf
        self.best_params = None
        self.last_sent = time.time()
