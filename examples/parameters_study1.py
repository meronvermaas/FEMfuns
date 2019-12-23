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

material = "dispersive" # options: "resistive" || "capacitive" || "dispersive"
frequencies = [1e7]

unit = 1e-2 #cm

'''
options:
	1) "name": [marker, conductivity]
	2) "name": [marker, admittivity]
	3) "name": [marker, colecolename] (in this case, name should be a material named in the cole cole table: help(disp.colecole_gabriel)
'''
volume_markers = {
'white': [32, 'Brain_Grey_Matter'],
'grey': [64, 'Brain_Grey_Matter'],
'csf': [96, 'Cerebro_Spinal_Fluid'],
'skull': [128, 'Bone_Cortical'],
'scalp': [160, 'Skin_Dry']
}

monopole_list = []
# add radial, tangential and 45-degree bipoles at increasing distances
for ii in np.arange(5):
    tmp_src_rad = [0., 0., 7.85- ii*.1, 1]
    tmp_snk_rad = [0., 0., 7.75- ii*.1, -1]
    tmp_src_tan = [0., -0.05, 7.8- ii*.1, 1]
    tmp_snk_tan = [0., 0.05, 7.8- ii*.1, -1]
    tmp_src_mix = [0., -0.0353, 7.835- ii*.1, 1]
    tmp_snk_mix = [0., 0.0353, 7.764- ii*.1, -1]

    monopole_list.append({'monopoles':[tmp_src_rad,tmp_snk_rad], 'name': 'rad'})
    monopole_list.append({'monopoles':[tmp_src_tan,tmp_snk_tan], 'name': 'tan'})
    monopole_list.append({'monopoles':[tmp_src_mix,tmp_snk_mix], 'name': 'mix'})

#########################################################################################################
####### coordinates as in https://github.com/Neuroinflab/fourspheremodel/blob/master/parameters.py ######
#########################################################################################################

theta, phi_angle = np.mgrid[0:180:1, -90:90:1]
theta = theta.flatten()
phi_angle = phi_angle.flatten()

theta_r = np.deg2rad(theta)
phi_angle_r = np.deg2rad(phi_angle)

rad_tol = 1e-2
brain_rad = 7.9
csftop_rad = 8.
skull_rad = 8.5
scalp_rad = 9.
x_points = (brain_rad - rad_tol) * np.sin(theta_r) * np.cos(phi_angle_r)
y_points = (brain_rad - rad_tol) * np.sin(theta_r) * np.sin(phi_angle_r)
z_points = (brain_rad - rad_tol) * np.cos(theta_r)

ele_coords = np.vstack((x_points, y_points, z_points)).T

