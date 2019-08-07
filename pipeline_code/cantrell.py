import numpy as np


def double_layer_impedance(frequency):
    #assuming doublelayer = 1e-9 #1nm

    if frequency == 0:
        frequency = 1e-16
        #return (338.,0.) #freq independent values from joucla
        #error('double layer impedance at frequency of zero can not be calculated, divide by 0')


    K = 1.57 #Ohm m2 s-beta
    beta = 0.91
    omega = 2*np.pi*frequency

    Zcpa = K*(1j*omega)**-beta # equation (3) in cantrell

    I0 = 6.41e-4 #Am-2
    T = 298 #K
    R = 8.3144598 #J mol-1 K-1
    F = 96485.33289 #C mol-1
    n = 2 #electrons/molecule

    Rct = R*T/(n*F*I0) # equation (5) in cantrell

    impedivity = Zcpa*Rct/(Zcpa+Rct)
    admittivity = 1./impedivity
    return np.real(admittivity), np.imag(admittivity)
