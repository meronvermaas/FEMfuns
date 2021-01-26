import numpy as np

'''
Volume conduction parameters
 - type of material (resistive, capacitive, dispersive)
 - type of electrode
 - monopole(s) list
 - marker link to material
 - measure coordinates
 - si unit
'''

material = "resistive" # "resistive" || "capacitive" || "dispersive"

unit = 1e-3 #mm

'''
options:
	1) "name": [marker, conductivity]
	2) "name": [marker, admittivity]
	3) "name": [marker, colecolename] (in this case, name should be a material named in the cole cole table: help(disp.colecole_gabriel)
'''
volume_markers = {
'white': [1, 0.25],
'grey': [2, 0.28],
'csf': [3, 1.59],
'skull': [4, 3.5e-3],
'scalp': [5, 0.17]
}
# add electrode region, purely insulating, use boundary_markers to describe surface conductance
for ii in range(6,70):
    volume_markers['elec_'+str(ii-5)] = [ii , 0]

volume_markers['grid'] = [70, 0]

'''
options:
	1) "name": [marker]
	2) "name": [marker, surface_conductance, 'ext'/'int']
		(surface_conductance can be capacitive or resistive)
		'ext'/'int' indicates interior or exterior boundary condition
	3) "name": [marker, potential, surface_conductance, 'ext'/'int']
'''
boundary_markers = {
#'white': [6],
'grey': [1]
#'csf': [8],
#'skull': [9] #,
#'scalp': [10]
}

for ii in range(2,66):
    boundary_markers['elec_'+str(ii-1)] = [ii ] #, 1e3, 'int']


monopole_list = []
#frequencies = [1000]
#coords_list = []
