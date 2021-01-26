import numpy as np
import matplotlib.pyplot as plt

max_all = []

nr_elecs = []
for ii in range(0,16):
    if ii == 0:
        lf_ds = np.load('results/patchb_4pointelec_M1S1_LFM_DS.npy')
    elif ii == 1:
        lf_ds = np.load('results/patchb_16pointelec_M1S1_LFM_DS.npy')
    elif ii == 2:
        lf_ds = np.load('results/patchb_64pointelec_M1S1_LFM_DS.npy')
    else:
        lf_ds = np.load('results/patchb_hd'+str(ii)+'pointelec_M1S1_LFM_DS.npy')
    max_all.append(np.amax(lf_ds,0))
    SVD_grid_DS = np.zeros([lf_ds.shape[1],lf_ds.shape[1]])

    for idx0 in range(lf_ds.shape[1]):
        for idx1 in range(lf_ds.shape[1]):
            U,S,Vt = np.linalg.svd([lf_ds[:,idx0],lf_ds[:,idx1]],full_matrices=False)
            if np.amax(lf_ds[:,idx0]) > 1 and np.amax(lf_ds[:,idx1]) > 1:
                SVD_grid_DS[idx0,idx1] = S[1]/S[0]
            else:
                SVD_grid_DS[idx0,idx1] = np.nan
    #deep x superficial
    SVD_grid_DS1 = SVD_grid_DS[1::2,::2]
    #nr of electrodes
    nr_elecs.append(lf_ds.shape[0])

nr_elecs = np.asarray(nr_elecs)

plt.rcParams.update({'font.size': 12})
n = 10
colors1 = plt.cm.gray(np.linspace(0., 1, n))
colors2 = plt.cm.spring(np.linspace(0, 1, n))
# combine them and build a new colormap
colors = np.vstack((colors1[:5,:], colors2[5:,:]))
labels = ['M1 d1','M1 d2','M1 d3','M1 d4','M1 d5','S1 d1','S1 d2','S1 d3','S1 d4','S1 d5']

##########plot everything
max_all = np.asarray(max_all)
max_all_sorted = max_all[np.argsort(nr_elecs),:]

odd = max_all_sorted[:,1::2]
even = max_all_sorted[:,::2]
odd[:,:5] = np.fliplr(odd[:,:5])
odd[:,5:] = np.fliplr(odd[:,5:])
max_all_sorted[:,1::2] = odd
even[:,:5] = np.fliplr(even[:,:5])
even[:,5:] = np.fliplr(even[:,5:])
max_all_sorted[:,::2] = even

plt.figure()
clridx = 0
for y_arr, label in zip(max_all_sorted[:,1::2].T, labels):
    plt.plot(nr_elecs[np.argsort(nr_elecs)], y_arr, 'o-', label=label, color=colors[clridx])
    clridx += 1
plt.xlabel('# of electrodes')
plt.ylabel('maximum $\mu$V')
plt.title('Deep patches')
plt.hlines(y = 1., xmin = nr_elecs.min(), xmax = nr_elecs.max(), color ='r',linewidth=3,linestyles ='dashed',label='threshold')
plt.xticks([nr_elecs.min()] + list(plt.xticks()[0]))
plt.autoscale(enable=True, axis='both', tight=True)

plt.figure()
clridx = 0
for y_arr, label in zip(max_all_sorted[:,::2].T, labels):
    plt.plot(nr_elecs[np.argsort(nr_elecs)], y_arr, '-o',label=label, color=colors[clridx])
    clridx += 1
plt.xlabel('# of electrodes')
plt.ylabel('maximum $\mu$V')
plt.title('Superficial patches')
plt.hlines(y = 1., xmin = nr_elecs.min(), xmax = nr_elecs.max(), color ='r',linewidth=3,linestyles ='dashed',label='threshold')
plt.xticks([nr_elecs.min()] + list(plt.xticks()[0]))
plt.autoscale(enable=True, axis='both', tight=True)

plt.show()

