import numpy as np


def double_layer_impedance(frequency):
    #assuming doublelayer = 1e-9 #1nm

    frequency = np.where(frequency==0,1e-20,frequency)

    I0 = 6.41e-4 #Am-2
    T = 298 #K
    R = 8.3144598 #J mol-1 K-1
    F = 96485.33289 #C mol-1
    n = 2 #electrons/molecule

    Rct = R*T/(n*F*I0) # equation (5) in cantrell

    K = 1.57 #Ohm m2 s-beta
    beta = 0.91
    omega = 2*np.pi*frequency

    Zcpa = K*(1j*omega)**-beta # equation (3) in cantrell

    impedivity = (Zcpa*Rct)/(Zcpa+Rct)
    admittivity = 1./impedivity

    return admittivity
