import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TKAgg')
from FEMfuns import *

dispcap = np.load('results/fem2_capdisp.npy')

# get the index of a vertex on the stimulating electrode surface
mesh = Mesh("mesh/joucla_correct_cells.xml")
boundaries = MeshFunction("size_t", mesh, "mesh/joucla_correct_cells_facet_region.xml")
V = FunctionSpace(mesh, "CG", 1)
gdim = mesh.geometry().dim()
dofs = V.tabulate_dof_coordinates().reshape((-1, gdim))
stim = 68
itstim = SubsetIterator(boundaries,stim)

for fs in itstim:
    continue
stim_idx = np.argmin(np.sum(np.abs(dofs-fs.midpoint()[:]),axis=1))

xlim = 4e-4
t = dispcap[0][1][4]
x = dispcap[0][1][5]

unit = 1e3 #s to ms
plt.rcParams.update({'font.size': 14})
plt.figure()
plt.title('Interface stimulating electrode vertex potential')
plt.plot(t*unit,x, label='pulse')
plt.plot(t*unit,dispcap[0][0][:,stim_idx].real,label='dispersive tissue, electrode RC',linewidth=4)
plt.plot(t*unit,dispcap[1][0][:,stim_idx].real,label='dispersive tissue, electrode pseudocapacitive',linewidth=3)
plt.plot(t*unit,dispcap[2][0][:,stim_idx].real,label='capacitive tissue, electrode RC',linewidth=2)
plt.plot(t*unit,dispcap[3][0][:,stim_idx].real,label='capacitive tissue, electrode pseudocapacitive')
plt.legend(loc=3, prop={'size': 13})
plt.xlim([-xlim*unit,xlim*unit])
plt.xlabel('Time (ms)')
plt.ylabel('Volts')
plt.savefig('fig2_dispcap.pdf')
