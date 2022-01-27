from FEMfuns import *
from os.path import exists
import sys
import params as params
from scipy.io import savemat

# total arguments
n = len(sys.argv)

# check that xml mesh is given and domains exist
if (n < 3 and not exists(sys.argv[1][0:-4]+'_physical_region.xml')) or not exists(sys.argv[1]):
    raise ValueError('Mesh filename or subdomains file could not be found')

mesh_filename = sys.argv[1]
if n==3 and sys.argv[2][-20:] == '_physical_region.xml':
    mesh_materials_filename = sys.argv[2]
else:
    mesh_materials_filename = sys.argv[1][0:-4]+'_physical_region.xml'
if n==4 and sys.argv[3][-17:] == '_facet_region.xml':
    mesh_boundaries_filename = sys.argv[3]
elif exists(sys.argv[1][0:-4]+'_facet_region.xml'):
    mesh_boundaries_filename = sys.argv[1][0:-4]+'_facet_region.xml'
else:
    mesh_boundaries_filename = None

FEMsim = FEM_simulation(params, mesh_filename=mesh_filename, mesh_materials_filename=mesh_materials_filename,mesh_boundaries_filename=mesh_boundaries_filename)

FEMsim.main(solver_method='cg', preconditioner='ilu')

surfintmarker = np.array([], dtype=int)
for ii in params.boundary_markers.values():
    if len(ii) == 1:
        surfintmarker = np.append(surfintmarker, ii)
    else:
        surfintmarker = np.append(surfintmarker, ii[0])

if len(surfintmarker)>0:
    elvals = FEMsim.get_poisson_values(surfintmarker,allsols=True)
else:
    elvals = FEMsim.get_poisson_values(params.elecpos,allsols=True)


savemat("femfuns_lfmatrix.mat", {'lf':elvals})
#print(elvals)

if 'brain' in params.volume_markers.keys():
    FEMsim.save_pvds(filename = 'brain_solution_example.pvd', subdomain_marker=params.volume_markers['brain'][0]) #allsols=True optional to get export all source solutions
elif 'grey' in params.volume_markers.keys():
    FEMsim.save_pvds(filename = 'grey_solution_example.pvd', subdomain_marker=params.volume_markers['grey'][0])
else:
    FEMsim.save_pvds(filename = 'solution_example.pvd')
