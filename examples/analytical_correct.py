#######################################################################################################
######taken from https://github.com/Neuroinflab/fourspheremodel/blob/master/analytical_correct.py######
######adjustments: alternative parameters, 1000 terms series expansion, calculate at brain surface#####
#######################################################################################################

import os
import numpy as np
from scipy.special import lpmv
import parameters_study1 as params

import argparse


def V(n):
    k = (n+1.) / n
    Factor = ( ( r34**n - (r43**(n+1)) ) / ( (k*(r34**n)) + (r43**(n+1)) ) )
    num = (s34/k) - Factor
    den = s34 + Factor
    return (num / den)


def Y(n):
    k = n / (n+1.)
    Factor = ( ( (r23**n) * k) - V(n)*(r32**(n+1))) / (r23**n + V(n)*(r32**(n+1)))
    num = (s23*k) - Factor
    den = s23 + Factor
    return (num / den)


def Z(n):
    k = (n+1.) / n
    num = (r12**n - k*Y(n)*(r21**(n+1)) ) / (r12**n + Y(n)*(r21**(n+1)))
    return num


def A1(n):
    num = (rz1**(n+1))* (Z(n) + s12*((n+1.)/n))
    den = s12 - Z(n)
    return num / den


def A2(n):
    num = A1(n) + (rz1**(n+1))
    den = (Y(n)*(r21**(n+1))) + r12**n
    return num / den


def B2(n):
    return A2(n)*Y(n)


def A3(n):
    num = A2(n) + B2(n)
    den = r23**n + (V(n)*(r32**(n+1)))
    return num / den


def B3(n):
    return A3(n)*V(n)


def A4(n):
    num = A3(n) + B3(n)
    k = (n+1.) / n
    den = (k*(r34**n)) + (r43**(n+1))
    return k*(num / den)


def B4(n):
    return A4(n)* (n / (n+1.))


def H(n, r_ele=params.brain_rad):
    if r_ele < params.brain_rad:
        T1 = ((r_ele / params.brain_rad)**n) * A1(n)
        T2 = ((rz / r_ele)**(n + 1))
    elif r_ele < params.csftop_rad:
        T1 = ((r_ele / params.csftop_rad)**n) * A2(n)
        T2 = ((params.csftop_rad / r_ele)**(n + 1)) * B2(n)
    elif r_ele < params.skull_rad:
        T1 = ((r_ele / params.skull_rad)**n) * A3(n)
        T2 = ((params.skull_rad / r_ele)**(n + 1)) * B3(n)
    elif r_ele <= params.scalp_rad:
        T1 = ((r_ele / params.scalp_rad)**n) * A4(n)
        T2 = ((params.scalp_rad / r_ele)**(n + 1)) * B4(n)
    else:
        print("Invalid electrode position")
        return
    return (T1 + T2)


def adjust_theta():
    ele_pos = params.ele_coords
    dp_loc = (np.array(src_pos) + np.array(snk_pos)) / 2.
    ele_dist = np.linalg.norm(ele_pos, axis=1)
    dist_dp = np.linalg.norm(dp_loc)
    cos_theta = np.dot(ele_pos, dp_loc) / (ele_dist * dist_dp)
    cos_theta = np.nan_to_num(cos_theta)
    theta = np.arccos(cos_theta)
    return theta


def adjust_phi_angle(p):
    ele_pos = params.ele_coords
    r_ele = np.sqrt(np.sum(ele_pos ** 2, axis=1))
    dp_loc = (np.array(src_pos) + np.array(snk_pos)) / 2.
    proj_rxyz_rz = (np.dot(ele_pos, dp_loc) / np.sum(dp_loc **2)).reshape(len(ele_pos),1) * dp_loc.reshape(1, 3)
    rxy = ele_pos - proj_rxyz_rz
    x = np.cross(p, dp_loc)
    cos_phi = np.dot(rxy, x.T) / np.dot(np.linalg.norm(rxy, axis=1).reshape(len(rxy),1), np.linalg.norm(x, axis=1).reshape(1, len(x)))
    cos_phi = np.nan_to_num(cos_phi)
    phi_temp = np.arccos(cos_phi)
    phi = phi_temp
    range_test = np.dot(rxy, p.T)
    for i in range(len(r_ele)):
        for j in range(len(p)):
            if range_test[i, j] < 0:
                phi[i,j] = 2 * np.pi - phi_temp[i, j]
    return phi.reshape(180 * 180)


def decompose_dipole(I):
    P = np.array([np.array(src_pos) * I - np.array(snk_pos) * I])
    dp_loc = (np.array(src_pos) + np.array(snk_pos)) / 2.
    dist_dp = np.linalg.norm(dp_loc)
    dp_rad = (np.dot(P, dp_loc) / dist_dp) * (dp_loc / dist_dp)
    dp_tan = P - dp_rad
    return P, dp_rad, dp_tan


def conductivity(sigma_skull):
    s12 = params.sigma_brain / params.sigma_csf
    s23 = params.sigma_csf / sigma_skull
    s34 = sigma_skull / params.sigma_scalp
    return s12, s23, s34


def compute_phi(s12, s23, s34, I):
    P, dp_rad, dp_tan = decompose_dipole(I)
    adjusted_theta = adjust_theta()

    adjusted_phi_angle = adjust_phi_angle(dp_tan)  # params.phi_angle_r

    dp_loc = (np.array(src_pos) + np.array(snk_pos)) / 2
    sign_rad = np.sign(np.dot(P, dp_loc))
    mag_rad = sign_rad * np.linalg.norm(dp_rad)
    mag_tan = np.linalg.norm(dp_tan)  # sign_tan * np.linalg.norm(dp_tan)

    coef = H(n)
    cos_theta = np.cos(adjusted_theta)

    # radial
    n_coef = n * coef
    rad_coef = np.insert(n_coef, 0, 0)
    Lprod = np.polynomial.legendre.Legendre(rad_coef)
    Lfactor_rad = Lprod(cos_theta)
    rad_phi = mag_rad * Lfactor_rad

    # #tangential
    Lfuncprod = []
    for tt in range(params.theta_r.size):
        Lfuncprod.append(np.sum([C * lpmv(1, P_val, cos_theta[tt])
                                 for C, P_val in zip(coef, n)]))
    tan_phi = -1 * mag_tan * np.sin(adjusted_phi_angle) * np.array(Lfuncprod)

    return (rad_phi + tan_phi) / (4 * np.pi * params.sigma_brain * (rz**2))


parser = argparse.ArgumentParser()
parser.add_argument('--directory', '-d',
                    default='results',
                    dest='results',
                    help='a path to the result directory')

args = parser.parse_args()

if not os.path.exists(args.results):
    os.makedirs(args.results)

# scalp_rad = scalp_rad - rad_tol
#rz = params.dipole_loc
#rz1 = rz / params.brain_rad
#r12 = params.brain_rad / params.csftop_rad
#r23 = params.csftop_rad / params.skull_rad
#r34 = params.skull_rad / params.scalp_rad

#r1z = 1. / rz1
#r21 = 1. / r12
#r32 = 1. / r23
#r43 = 1. / r34

I = 1.
n = np.arange(1, 1000)

#####################################################################################################################
######main changes compared to https://github.com/Neuroinflab/fourspheremodel/blob/master/analytical_correct.py######
######occur below, to correctly set parameters#######################################################################
#####################################################################################################################

import dispersion_spectra as disp

omega = 2.*np.pi*params.frequencies[0]
e0 = 8.854187817e-12 #F/m

#brain matter
tmpeps, tmpsig = disp.colecole_gabriel(params.frequencies[0], 'Brain_Grey_Matter')
params.sigma_brain = tmpsig * params.unit + 1j * tmpeps * omega * e0 * params.unit

#csf
tmpeps, tmpsig = disp.colecole_gabriel(params.frequencies[0], 'Cerebro_Spinal_Fluid')
params.sigma_csf = tmpsig * params.unit + 1j * tmpeps * omega * e0 * params.unit

#skull
tmpeps, tmpsig = disp.colecole_gabriel(params.frequencies[0], 'Bone_Cortical')
params.sigma_skull = tmpsig * params.unit + 1j * tmpeps * omega * e0 * params.unit

#scalp
tmpeps, tmpsig = disp.colecole_gabriel(params.frequencies[0], 'Skin_Dry')
params.sigma_scalp = tmpsig * params.unit + 1j * tmpeps * omega * e0 * params.unit

allphi = []
for bipole in params.monopole_list:
    print('Now computing for dipole: ', bipole['name'])
    # changed from original, readout parameters with slightly different format
    src_pos = bipole['monopoles'][0][0:3]
    snk_pos = bipole['monopoles'][1][0:3]

    # added to original code to set for all dipole locations
    params.dipole_loc = np.round(bipole['monopoles'][0][2]-.005,1)
    rz = params.dipole_loc
    rz1 = rz / params.brain_rad
    r12 = params.brain_rad / params.csftop_rad
    r23 = params.csftop_rad / params.skull_rad
    r34 = params.skull_rad / params.scalp_rad

    r1z = 1. / rz1
    r21 = 1. / r12
    r32 = 1. / r23
    r43 = 1. / r34

    s12, s23, s34 = conductivity(params.sigma_skull)
    phi = compute_phi(s12, s23, s34, I)

    # changed from original code, commented out
    #s12, s23, s34 = conductivity(params.sigma_skull40)
    #phi_40 = compute_phi(s12, s23, s34, I)

    #s12, s23, s34 = conductivity(params.sigma_skull80)
    #phi_80 = compute_phi(s12, s23, s34, I)

    #s12 = s23 = s34 = 1.
    #phi_lim = compute_phi(s12, s23, s34, I)

    #with open(os.path.join(args.results,
    #                       'Analytical_' + dipole['name'] + '.npz'), 'wb') as f:
    #    np.savez(f, phi_20=phi_20, phi_40=phi_40, phi_80=phi_80, phi_lim=phi_lim)
    allphi.append(phi)

np.save(os.path.join(args.results,'Analytical_allphi'+str(int(params.frequencies[0]))+'.npy'), allphi)
