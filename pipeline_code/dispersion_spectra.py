import numpy as np
import os
this_dir, this_filename = os.path.split(__file__)

colecole_params = np.genfromtxt(this_dir+"/colecole_parameters.csv", delimiter=',', \
                  dtype = None, skip_header=1, \
                  names = ['Tissue','ef','del1','tau1_ps','alf1','del2','tau2_ns','alf2','sig','del3','tau3_us','alf3','del4','tau4_ms','alf4'])

def get_gabriel_params(tissue_type):
    '''
 Possible tissue types are:\n Aorta\n Bladder\n Blood\n Bone_Cancellous\n Bone_Cortical\n Bone_Marrow_Infiltrated\n Bone_Marrow_Not_Infiltrated\n Brain_Grey_Matter\n Brain_White_Matter\n Breast_fat\n Cartilage\n Cerebellum\n Cerebro_Spinal_Fluid\n Cervix\n Colon\n Cornea\n Dura\n Eye_Tissues_Sclera\n Fat_Average_Infiltrated\n Fat_Not_Infiltrated\n Gall_Bladder\n Gall_Bladder_Bile\n Heart\n Kidney\n Lens_Cortex\n Lens_Nucleus\n Liver\n Lung_Deflated\n Lung_Inflated\n Muscle\n Nerve\n Ovary\n Skin_Dry\n Skin_Wet\n Small_Intestine\n Spleen\n Stomach\n Tendon\n Testis\n Thyroid\n Tongue\n Trachea\n Uterus\n Vitreous_Humor
 returns epsiloninf, Depsilon[De1,De2,De3,De4], tau[t_ps,t_ns,t_us,t_ms], alpha[a1,a2,a3,a4],sigma
    '''

    tissue_arg = [i for i, tissues in enumerate(colecole_params) if tissue_type in tissues]
    if tissue_arg:
        vals = colecole_params[tissue_arg[0]]
        epsilon_inf = vals[1]
        Depsilon = [vals[2],vals[5],vals[9],vals[12]]
        tau = [vals[3]*1e-12,vals[6]*1e-9,vals[10]*1e-6,vals[13]*1e-3] #ps,ns,us,ms
        alpha = [vals[4],vals[7],vals[11],vals[14]]
        sigma = vals[8]

        return epsilon_inf, Depsilon, tau, alpha, sigma
    else:
        return "this tissue is not in the list"

def colecole_gabriel(frequency,tissue_type):
    '''
 Taken from Gabriel et al 1996, The dielectric properties of biological tissues: III. Parametric models for the dielectric spectrum of tissues
 Possible tissue types are:\n Aorta\n Bladder\n Blood\n Bone_Cancellous\n Bone_Cortical\n Bone_Marrow_Infiltrated\n Bone_Marrow_Not_Infiltrated\n Brain_Grey_Matter\n Brain_White_Matter\n Breast_fat\n Cartilage\n Cerebellum\n Cerebro_Spinal_Fluid\n Cervix\n Colon\n Cornea\n Dura\n Eye_Tissues_Sclera\n Fat_Average_Infiltrated\n Fat_Not_Infiltrated\n Gall_Bladder\n Gall_Bladder_Bile\n Heart\n Kidney\n Lens_Cortex\n Lens_Nucleus\n Liver\n Lung_Deflated\n Lung_Inflated\n Muscle\n Nerve\n Ovary\n Skin_Dry\n Skin_Wet\n Small_Intestine\n Spleen\n Stomach\n Tendon\n Testis\n Thyroid\n Tongue\n Trachea\n Uterus\n Vitreous_Humor
 returns epsilon, sigma
    '''

    if np.abs(frequency) < 100: #False:
        if tissue_type == 'Brain_Grey_Matter':
            return 0, .239
        elif tissue_type == 'Brain_White_Matter':
            return 0, .265
        elif tissue_type == 'Cerebro_Spinal_Fluid':
            return 0, 1.78
        elif tissue_type == 'Bone_Cortical':
            return 0, 3.5e-3
        elif tissue_type == 'Skin_Wet' or tissue_type == 'Skin_Dry':
            return 0, 3.2e-5
        elif tissue_type == 'Cerebellum':
            return 0, .66
        else:
            return "This material conductivity below 100Hz is unknown"
    else:
        e_0 = 8.85e-12 #F/m
        omega = frequency*2*np.pi

        # get params from gabriel
        epsilon_inf, Depsilon, tau, alpha, sigma = get_gabriel_params(tissue_type)

        # calc colecole
        eps = epsilon_inf \
              + Depsilon[0] / (1+(1j*omega*tau[0])**(1-alpha[0])) \
              + Depsilon[1] / (1+(1j*omega*tau[1])**(1-alpha[1])) \
              + Depsilon[2] / (1+(1j*omega*tau[2])**(1-alpha[2])) \
              + Depsilon[3] / (1+(1j*omega*tau[3])**(1-alpha[3])) \
              + sigma / (1j*omega*e_0)

        return np.real(eps), np.real(1j*e_0*omega*eps)

def dispersion_raicu_skin(frequency):
    '''
    Taken from Raicu et al 2000, A quantitative approach to the dielectric properties of the skin
    These should be more accurate than the properties Gabriel's approximation
    '''

    # params from raicu
    delta = 2200
    fc = 2.51e6 #Hz
    alpha = .21
    beta = .152
    gamma = 1
    kl = 32e-6 #S/m
    eh = 61

    e_0 = 8.85e-12 #F/m
    omega = frequency*2*np.pi

    eps =  eh \
           + delta / ( (1j*frequency/fc)**alpha + (1j*frequency/fc)**(1-beta)  )**gamma \
           + kl / (1j*omega*e_0)

    return  np.real(eps), np.real(1j*e_0*omega*eps)
