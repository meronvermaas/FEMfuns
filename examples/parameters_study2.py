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

unit = 1e-3 #mm

'''
options:
	1) "name": [marker, conductivity]
	2) "name": [marker, admittivity]
	3) "name": [marker, colecolename] (in this case, name should be a material named in the cole cole table: help(disp.colecole_gabriel)
'''
volume_markers = {
'brain': [1, 'Brain_Grey_Matter'],
'csf': [2, 'Cerebro_Spinal_Fluid']
}


'''
options:
	1) "name": [marker]
	2) "name": [marker, surface_conductance, 'ext'/'int']
		(surface_conductance can be capacitive or resistive, i.e. a real or complex number
                 OR the string 'dispersive' to indicate cantrell (2009) should be used)
		'ext'/'int' indicates interior or exterior boundary condition
	3) "name": [marker, potential, surface_conductance, 'ext'/'int']
'''
boundary_markers = {
'stim_elec': [68, 1, 'dispersive', 'ext'],
'ground_elec': [1, 0, 'dispersive', 'ext']
}
