from FEMfuns import *
import parameters_study2 as params

FEMsim = FEM_simulation(params, mesh_filename='mesh/joucla_correct_cells.xml', mesh_materials_filename='mesh/joucla_correct_cells_physical_region.xml', mesh_boundaries_filename='mesh/joucla_correct_cells_facet_region.xml')

FEMsim.make_pulse(1,2e-4,100001)

######################################################
###dispersive tissue, pseudocapacitive electrode######
######################################################
FEMsim.main(solver_method='gmres', preconditioner='hypre_euclid', relative_tolerance=1e-2)

######################################################
###dispersive tissue, capacitive electrode############
######################################################
# take parallel circuit admittance at average frequency
# real surface conductance in S/m2 and capacitance in F/m2 (will be multiplied by omega to get susceptance in Siemens)
# these values serve as an example, any complex number describing the RC of the interface per square meter suffices
surface_admittance = ct.double_layer_impedance(FEMsim.params.frequencies.max()/2)
surface_admittance = complex(surface_admittance.real, surface_admittance.imag/(np.pi*2*(FEMsim.params.frequencies.max()/2)))

for bound in FEMsim.params.boundary_markers:
    FEMsim.params.boundary_markers[bound][2] = surface_admittance
FEMsim.main(solver_method='gmres', preconditioner='hypre_euclid', relative_tolerance=1e-2)


######################################################
###capacitive tissue, capacitive electrode############
######################################################
FEMsim.params.material = 'capacitive'
brain_epsilon, brain_sigma = disp.colecole_gabriel(FEMsim.params.frequencies.max()/2, 'Brain_Grey_Matter')
brain_admittivity = complex(brain_sigma,brain_epsilon)
csf_eps, csf_sig = disp.colecole_gabriel(FEMsim.params.frequencies.max()/2, 'Cerebro_Spinal_Fluid')
csf_admittivity = complex(csf_sig, csf_eps)

FEMsim.params.volume_markers['brain'][1] = brain_admittivity
FEMsim.params.volume_markers['csf'][1] = csf_admittivity
FEMsim.geometry.init_admittivity(FEMsim)

FEMsim.main(solver_method='gmres', preconditioner='hypre_euclid', relative_tolerance=1e-2)

############################################################
###capacitive tissue, pseudocapacitive electrode############
############################################################

for bound in FEMsim.params.boundary_markers:
    FEMsim.params.boundary_markers[bound][2] = 'dispersive'
FEMsim.main(solver_method='gmres', preconditioner='hypre_euclid', relative_tolerance=1e-2)

#plot a vertex to inspect
#for ii in range(len(FEMsim.solutions)):
#    plt.plot(FEMsim.solutions[ii][1][4], FEMsim.solutions[ii][0][:,1000].real)
#plt.show()