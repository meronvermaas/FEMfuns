from FEMfuns import *
import parameters_pointelec as params

FEMsim = FEM_simulation(params, mesh_filename='mesh/headgeom.xml', mesh_materials_filename='mesh/headgeom_physical_region.xml', mesh_boundaries_filename='mesh/headgeom_facet_region.xml')

elec_coords = np.loadtxt('mesh/elec-coor-innerskull.txt')
M1_1 = np.loadtxt('sources/M1_finger1_dipole_coords_patches2.txt', delimiter=",")
M1_2 = np.loadtxt('sources/M1_finger2_dipole_coords_patches2.txt', delimiter=",")
M1_3 = np.loadtxt('sources/M1_finger3_dipole_coords_patches2.txt', delimiter=",")
M1_4 = np.loadtxt('sources/M1_finger4_dipole_coords_patches2.txt', delimiter=",")
M1_5 = np.loadtxt('sources/M1_finger5_dipole_coords_patches2.txt', delimiter=",")
S1_1 = np.loadtxt('sources/S1_finger1_dipole_coords_patches2.txt', delimiter=",")
S1_2 = np.loadtxt('sources/S1_finger2_dipole_coords_patches2.txt', delimiter=",")
S1_3 = np.loadtxt('sources/S1_finger3_dipole_coords_patches2.txt', delimiter=",")
S1_4 = np.loadtxt('sources/S1_finger4_dipole_coords_patches2.txt', delimiter=",")
S1_5 = np.loadtxt('sources/S1_finger5_dipole_coords_patches2.txt', delimiter=",")

grey_surf_marker = FEMsim.params.boundary_markers['grey'][0]

#M1 plane
plane_p1 = np.array([37.4,-3.8,50.8])
plane_n1 = np.array([0,0,1])
#S1 plane
plane_p2 = np.array([40.5,-11.4,44.8])
plane_n2 = np.array([0.577733778841725,0.281575490316184,0.766119392809281])

patchname = ['M1_1','M1_2','M1_3','M1_4','M1_5','S1_1','S1_2','S1_3','S1_4','S1_5']
patchidx = 0
for finger_dipoles in [M1_1,M1_2,M1_3,M1_4,M1_5,S1_1,S1_2,S1_3,S1_4,S1_5]:
    FEMsim.get_closest_boundary(grey_surf_marker, finger_dipoles)
    monopoles_deep = []
    monopoles_shallow = []
    dipole_strength = 1
    for idx,coord in enumerate(finger_dipoles):
        src_coord = coord+FEMsim.closest_normals[idx][:]*.5
        snk_coord = coord-FEMsim.closest_normals[idx][:]*.5
        if patchidx <= 4:
            planecheck = np.dot(plane_n1, coord-plane_p1)
        else:
            planecheck = np.dot(plane_n2, coord-plane_p2)

        if np.sign(planecheck) == 1:
            monopoles_shallow.append(np.append(src_coord,dipole_strength))
            monopoles_shallow.append(np.append(snk_coord,-dipole_strength))
        else:
            monopoles_deep.append(np.append(src_coord,dipole_strength))
            monopoles_deep.append(np.append(snk_coord,-dipole_strength))
    monopoles_shallow = np.asarray(monopoles_shallow)
    monopoles_deep = np.asarray(monopoles_deep)
    monopoles_shallow[:,-1] = monopoles_shallow[:,-1]/monopoles_shallow.shape[0]
    monopoles_deep[:,-1] = monopoles_deep[:,-1]/monopoles_deep.shape[0]
    FEMsim.params.monopole_list.append({'monopoles': monopoles_shallow, 'name':patchname[patchidx]+'_deep'})
    FEMsim.params.monopole_list.append({'monopoles': monopoles_deep, 'name':patchname[patchidx]+'_shallow'})
    patchidx += 1

tic = time()
FEMsim.main(solver_method='cg', preconditioner='ilu')
print('total calculation of all dipoles took ', (time()-tic)/60, ' minutes')

elvals_point = FEMsim.get_poisson_values(elec_coords[[6,14,16,46],:],allsols=True)
elvals_point = np.asarray(elvals_point)
###referencing
elvals_point = elvals_point.T - np.mean(elvals_point,1)
np.save('results/patchb_4pointelec_M1S1_LFM_DS.npy',elvals_point)
np.savetxt('results/patchb_4pointelec_M1S1_LFM_DS..txt',np.concatenate((elec_coords[[6,14,16,46],:],elvals_point),axis=1),delimiter=",")

elvals_point = FEMsim.get_poisson_values(elec_coords[[0,6,11,14,16,19,24,26,30,34,38,43,44,46,53,62],:],allsols=True)
elvals_point = np.asarray(elvals_point)
###referencing
elvals_point = elvals_point.T - np.mean(elvals_point,1)
np.save('results/patchb_16pointelec_M1S1_LFM_DS.npy',elvals_point)
np.savetxt('results/patchb_16pointelec_M1S1_LFM_DS.txt',np.concatenate((elec_coords[[0,6,11,14,16,19,24,26,30,34,38,43,44,46,53,62],:],elvals_point),axis=1),delimiter=",")

elvals_point = FEMsim.get_poisson_values(elec_coords,allsols=True)
elvals_point = np.asarray(elvals_point)
###referencing
elvals_point = elvals_point.T - np.mean(elvals_point,1)
np.save('results/patchb_64pointelec_M1S1_LFM_DS.npy',elvals_point)
np.savetxt('results/patchb_64pointelec_M1S1_LFM_DS.txt',np.concatenate((elec_coords,elvals_point),axis=1),delimiter=",")


for ii in range(1,16):
    elec_coords = np.loadtxt('mesh/HHD'+str(ii)+'_elec-coor-innerskull.txt', delimiter=",")
    elvals_point = FEMsim.get_poisson_values(elec_coords,allsols=True)
    elvals_point = np.asarray(elvals_point)
    ###referencing
    elvals_point = elvals_point.T - np.mean(elvals_point,1)
    np.save('results/patchb_hd'+str(ii)+'pointelec_M1S1_LFM_DS.npy',elvals_point)
    np.savetxt('results/patchb_hd'+str(ii)+'_pointelec_M1S1_LFM_DS.txt',np.concatenate((elec_coords,elvals_point),axis=1),delimiter=",")
