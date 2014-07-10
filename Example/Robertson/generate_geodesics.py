"""
Calculate geodesics for our example model
"""
from SloppyCell.ReactionNetworks import *
from numpy import *
import geodesics
import example_model

x = log(example_model.popt)
v = array([0,1])
tmax = 10
Avv = None

func = lambda logp: array(example_model.m.res_log_params(logp))
jacobian = lambda logp: array(example_model.m.jacobian_log_params_sens(logp).values())

u, vects = Utility.eig(example_model.jtj)
# Draw geodesics on equally spaced angles coming out from optimal parameter
# values, starting at sloppy direction.
theta0 = arctan2(vects[1,-1], vects[0,-1])
theta_list = theta0 + linspace(0,2*pi,8,endpoint=False)

result_list = []
for theta in theta_list:
    v = array([cos(theta), sin(theta)])
    xs, vs, ts = geodesics.geodesic(x, v, tmax, func, jacobian, Avv, 
                                    maxsteps=int(2e4), rtol=1e-3, atol=1e-3)
    result_list.append((xs,vs,ts))
    print('Got to t={0}'.format(ts[-1]))

Utility.save(result_list, 'example.geodesics.bpkl')
