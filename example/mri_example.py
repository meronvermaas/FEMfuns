import sys
appfolder = '/path/to/forward_model/folder')
sys.path.insert(0, appfolder)

from fw_fenics_classes import *
import parameters_mri_example as params

FEMsim = FEM_simulation(params, mesh_filename='mesh/almi_2elechole.xml', mesh_materials_filename='mesh/almi_2elechole_physical_region.xml', mesh_boundaries_filename='mesh/almi_2elechole_facet_region.xml')

#create one radial source at 1cm from elec somewhere in the grey matter
# more info help(FEM_simulation.source_locations)
FEMsim.source_locations(FEMsim.params.volume_markers['grey'][0],\
                        FEMsim.geometry.boundaries, \
                        FEMsim.params.boundary_markers['elec1'][0], \
                        1., 1., 1, ['rad'])
FEMsim.main()

FEMsim.set_dispersive_params()

FEMsim.main(solver_method='gmres', preconditioner='ilu')

FEMsim.save_pvds(allsols=True)
