import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TKAgg')

all_fem_comp = np.load('results/fem_foursphere10000000.npy')
allphi = np.load('results/Analytical_allphi10000000.npy')
diff = allphi-all_fem_comp
MRD = np.mean(np.abs(diff)/np.abs(np.amax(allphi,axis=1)[:,None]), axis=1) # mean relative difference
MRD.resize([5,3])

dists = np.arange(1,6) # closest distance to grey surf in mm

plt.rcParams.update({'font.size': 14})
plt.plot(dists, MRD[:,0],'-o', label='Radial dipole', linewidth=4, c='0.', markersize=8)
plt.plot(dists, MRD[:,1],'-v', label='Tangential dipole', linewidth=4, c='.4', markersize=8)
plt.plot(dists, MRD[:,2],'-d', label='45 dipole', linewidth=4, c='.8', markersize=8)
plt.legend(loc=1, prop={'size': 13})
plt.xlabel('Smallest distance to brain surface (mm)')
plt.ylabel('Mean relative difference')
plt.savefig('fig1_foursphere.pdf')
