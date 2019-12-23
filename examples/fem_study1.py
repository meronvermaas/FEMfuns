from FEMfuns import *
import parameters_study1 as params

FEMsim = FEM_simulation(params, mesh_filename='../mesh/sphere_4.xml', mesh_materials_filename='../mesh/sphere_4_physical_region.xml', mesh_boundaries_filename='../mesh/sphere_4_facet_region.xml', order=2)

FEMsim.main(solver_method='gmres', preconditioner='ilu')

fem_vals = FEMsim.get_poisson_values(FEMsim.params.ele_coords, allsols=True)

np.save('results/fem_foursphere'+str(int(params.frequencies[0]))+'.npy',fem_vals)
