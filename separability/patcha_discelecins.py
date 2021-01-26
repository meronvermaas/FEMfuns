from FEMfuns import *
import parameters_discelecins as params

FEMsim = FEM_simulation(params, mesh_filename='mesh/headgeom.xml', mesh_materials_filename='mesh/headgeom_physical_region.xml', mesh_boundaries_filename='mesh/headgeom_facet_region.xml')

elec_coords = np.loadtxt('mesh/elec-coor-innerskull.txt')
M1_1 = np.loadtxt('sources/M1_finger1_dipole_coords.txt', delimiter=",")
M1_2 = np.loadtxt('sources/M1_finger2_dipole_coords.txt', delimiter=",")
M1_3 = np.loadtxt('sources/M1_finger3_dipole_coords.txt', delimiter=",")
M1_4 = np.loadtxt('sources/M1_finger4_dipole_coords.txt', delimiter=",")
M1_5 = np.loadtxt('sources/M1_finger5_dipole_coords.txt', delimiter=",")
S1_1 = np.loadtxt('sources/S1_finger1_dipole_coords.txt', delimiter=",")
S1_2 = np.loadtxt('sources/S1_finger2_dipole_coords.txt', delimiter=",")
S1_3 = np.loadtxt('sources/S1_finger3_dipole_coords.txt', delimiter=",")
S1_4 = np.loadtxt('sources/S1_finger4_dipole_coords.txt', delimiter=",")
S1_5 = np.loadtxt('sources/S1_finger5_dipole_coords.txt', delimiter=",")

grey_surf_marker = FEMsim.params.boundary_markers['grey'][0]

patchname = ['M1_1','M1_2','M1_3','M1_4','M1_5','S1_1','S1_2','S1_3','S1_4','S1_5']
patchidx = 0
for finger_dipoles in [M1_1,M1_2,M1_3,M1_4,M1_5,S1_1,S1_2,S1_3,S1_4,S1_5]:
    FEMsim.get_closest_boundary(grey_surf_marker, finger_dipoles)
    monopoles = []
    dipole_strength = 1/finger_dipoles.shape[0]
    for idx,coord in enumerate(finger_dipoles):
        src_coord = coord+FEMsim.closest_normals[idx][:]*.5
        snk_coord = coord-FEMsim.closest_normals[idx][:]*.5
        monopoles.append(np.append(src_coord,dipole_strength))
        monopoles.append(np.append(snk_coord,-dipole_strength))
    FEMsim.params.monopole_list.append({'monopoles': monopoles, 'name':patchname[patchidx]})
    patchidx += 1

tic = time()
FEMsim.main(solver_method='cg', preconditioner='ilu')
print('total calculation of all dipoles took ', (time()-tic)/60, ' minutes')

elvals_disc = FEMsim.get_poisson_values(np.arange(2,66),allsols=True)
elvals_disc = np.asarray(elvals_disc)

#referencing
elvals_disc = elvals_disc.T - np.mean(elvals_disc,1)

np.save('results/patcha_discelecins_M1S1_LFM.npy',elvals_disc)
np.savetxt('results/patcha_discelecins_M1S1_LFM.txt',np.concatenate((elec_coords,elvals_disc),axis=1),delimiter=",")
