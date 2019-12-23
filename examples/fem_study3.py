from FEMfuns import *
import parameters_study3 as params

FEMsim = FEM_simulation(params, mesh_filename='mesh/ernie.xml', mesh_materials_filename='mesh/ernie_physical_region.xml', mesh_boundaries_filename='mesh/ernie_facet_region.xml')

FEMsim.main(solver_method='gmres', preconditioner='ilu', maximum_iterations=500)

FEMsim.save_pvds(allsols = True, filename = 'ernie_grey',subdomain_marker=2)
