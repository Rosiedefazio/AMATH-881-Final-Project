# import os

# import numpy as np
# import matplotlib.pyplot as plt


# x = np.linspace(0.0, 12.0, 600)
# y = np.exp(-0.12 * x) * np.sin(2.4 * x)

# fig, ax = plt.subplots(figsize=(7, 4.5))
# ax.plot(x, y, color='#1f77b4', linewidth=2.2)
# ax.set_title('Dummy Model A')
# ax.set_xlabel('time')
# ax.set_ylabel('response')
# ax.grid(alpha=0.3)

# savefig_path = os.getenv('PYDSTOOL_SAVEFIG')
# if savefig_path:
#     plt.savefig(savefig_path, dpi=160, bbox_inches='tight')
#     print('Saved plot to %s' % savefig_path)
# else:
#     plt.show()


""" EXAMPLE: Logistic map

    Drew LaMar, March 2006
"""

import os

from PyDSTool import *

pars = {'r': 0.}

icdict = {'x': 0.}

# Set up model
xstr = 'r*x*(1-x)'

DSargs = args(name='LogisticMap')
DSargs.pars = pars
DSargs.varspecs = {'x': xstr}
DSargs.ics = icdict
DSargs.ttype = int
DSargs.pdomain = {'r': [0.0, 4.0]}

testDS = Generator.MapSystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)


def orient_branch(point, label, parname='r', increasing=True):
    """Normalize branch direction so continuation follows a stable parameter orientation.

    The tangent returned by PyCont can flip sign depending on the linear algebra
    backend, so orient it explicitly along the free parameter used in this example.
    """

    branch = point.labels[label]['data'].branch
    want_positive = 1 if increasing else -1
    if branch[parname] * want_positive < 0:
        return {k: -v for k, v in branch.items()}
    return branch


PCargs = args(name='FP1', type='FP-C')
PCargs.freepars = ['r']
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 70
PCargs.MaxStepSize = 1e-1
PCargs.verbosity = 2
PCargs.LocBifPoints = 'all'
PCargs.StopAtPoints = 'B'
PCargs.SaveEigen = True
PCargs.SPOut = {'r': [0.1919191, 1.5353535]}
PyCont.newCurve(PCargs)

print('Computing curve...')
start = perf_counter()
PyCont['FP1'].forward()
print('done in %.3f seconds!' % (perf_counter()-start))

PCargs.name = 'FP2'
PCargs.initpoint = 'FP1:BP1'
PCargs.initdirec = orient_branch(PyCont['FP1'].getSpecialPoint('BP1'), 'BP')
PCargs.MaxNumPoints=50
PCargs.LocBifPoints = ['PD', 'B']
PyCont.newCurve(PCargs)

print('Computing second branch...')
start = perf_counter()
PyCont['FP2'].forward()
print('done in %.3f seconds!' % (perf_counter()-start))

PCargs.name = 'FP3'
PCargs.initpoint = 'FP2:PD1'
PCargs.initdirec = orient_branch(PyCont['FP2'].getSpecialPoint('PD1'), 'PD')
PCargs.MaxNumPoints = 40
PCargs.LocBifPoints = ['PD', 'B']
PCargs.period = 2
PyCont.newCurve(PCargs)

print('Computing 2-cycle branch...')
start = perf_counter()
PyCont['FP3'].forward()
PyCont['FP3'].backward()
PyCont['FP3'].cleanLabels()
print('done in %.3f seconds!' % (perf_counter()-start))

PCargs.name='FP4'
PCargs.initpoint = 'FP3:PD1'
PCargs.initdirec = orient_branch(PyCont['FP3'].getSpecialPoint('PD1'), 'PD')
PCargs.period = 4
PyCont.newCurve(PCargs)

print('Computing 1st 4-cycle branch...')
start = perf_counter()
PyCont['FP4'].forward()
PyCont['FP4'].backward()
PyCont['FP4'].cleanLabels()
print('done in %.3f seconds!' % (perf_counter()-start))

PCargs.name = 'FP5'
PCargs.initpoint = 'FP3:PD2'
PCargs.initdirec = orient_branch(PyCont['FP3'].getSpecialPoint('PD2'), 'PD')
PyCont.newCurve(PCargs)

print('Computing 2nd 4-cycle branch...')
start = perf_counter()
PyCont['FP5'].forward()
PyCont['FP5'].backward()
PyCont['FP5'].cleanLabels()
print('done in %.3f seconds!' % (perf_counter()-start))

# Plot
PyCont.display(stability=True)
plt.xlim([1, 4])
PyCont.plot.toggleAll('off')
plt.title('Logistic map')

savefig_path = os.getenv('PYDSTOOL_SAVEFIG')
if savefig_path:
    plt.savefig(savefig_path, dpi=160, bbox_inches='tight')
    print('Saved plot to %s' % savefig_path)
else:
    show()
