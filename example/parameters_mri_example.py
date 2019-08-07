import numpy as np



'''
Volume conduction parameters
 - si unit
 - type of material (resistive vs capacitive)
 - type of electrode
 - monopole(s) list
 - marker link to material
 - measure coordinates
'''

material = "resistive" # "resistive" || "capacitive" || "dispersive"
frequencies = [1000]

unit = 1e-2 #convert to cm

sigma_brain = 1./300. #+ 300.j

# two options:
# "name": [marker, conductivity]
# or
# "name": [marker, colecolename] (in this case, name should be a material named in the cole cole table: help(disp.colecole_gabriel)
volume_markers = {
'cerebellum': [13, sigma_brain],
'ventricles': [11, 5*sigma_brain],
'white': [9, sigma_brain],
'grey': [7, sigma_brain],
'csf': [5, 5*sigma_brain],
'skull': [3, sigma_brain / 20.],
'scalp': [1, sigma_brain]
}

boundary_markers = {
'elec1': [15],
'elec2': [16]
}
