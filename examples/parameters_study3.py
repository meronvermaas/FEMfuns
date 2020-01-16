import numpy as np

'''
Volume conduction parameters
 - type of material (resistive, capacitive, dispersive)
 - type of electrode (optional, in boundary_markers)
 - monopole(s) list (optional)
 - marker link to material (volume_markers)
 - measure coordinates (optional)
 - si unit
'''

material = "resistive" # options: "resistive" || "capacitive" || "dispersive"

unit = 1e-3 #mm

'''
options:
	1) "name": [marker, conductivity]
	2) "name": [marker, admittivity]
	3) "name": [marker, colecolename] (in this case, name should be a material named in the cole cole table: help(disp.colecole_gabriel)
        4) "0": array of size [nr_cells,9] describing the anisotropic tensor
'''
volume_markers = {
'0': np.load('mesh/aniso_tensor.npy')
}

monopole_list = []
# select point location in paraview
src_x,src_y,src_z = 18.1176, -24.7466, 96.1343
# look for normal of closest grey matter boundary
# performed separately with:
# grey_boundary = 1002
#FEMsim.get_closest_boundary(np.array([[src_x,src_y,src_z]]),grey_boundary)
#normal=FEMsim.closest_normals[0][:]
normal_x,normal_y,normal_z = 0.44825004, -0.48938227,  0.74804873

src = [src_x,src_y,src_z, 1]
#distance between source and sink is 1 mm
snk_rad = [src_x-normal_x,src_y-normal_y,src_z-normal_z, -1]
snk_tan = [src_x+normal_z,src_y-normal_y,src_z-normal_x, -1]
monopole_list.append({'monopoles':[src,snk_rad], 'name': 'rad'})
monopole_list.append({'monopoles':[src,snk_tan], 'name': 'tan'})
