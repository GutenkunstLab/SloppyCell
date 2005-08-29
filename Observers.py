import scipy

def print_hess_elements(**args):
    if args['event'] == 'hessian element':
        elem = args['element']
        ii = args['i']
        jj = args['j']

        print 'hessian element %i, %i: %g' % (ii, jj, elem)

class CostPrinter:
    def __init__(self, skip=1):
        self.skip = skip
        self.reset()

    def __call__(self, **args):
        if args['event'] == 'evaluation':
            cost = args['cost']
            params = args['params']

            self.lowest_cost = min(cost, self.lowest_cost)
            if self.ii % self.skip == 0:
                print 'call %i: cost: %g, best so far: %g' % (self.ii, cost, 
                                                              self.lowest_cost)
            self.ii += 1

    def reset(self):
        self.ii = 0
        self.lowest_cost = scipy.inf

def print_all_costs(**args):
    if args['event'] == 'cost':
        print args['cost']
