import numpy as np
import matplotlib.pyplot as plt

lf = np.load('results/patcha_pointelec_M1S1_LFM.npy')

SVD_grid = np.zeros([lf.shape[1],lf.shape[1]])
for idx0 in range(lf.shape[1]):
    for idx1 in range(lf.shape[1]):
        U,S,Vt = np.linalg.svd([lf[:,idx0],lf[:,idx1]],full_matrices=False)
        if np.amax(np.abs(lf[:,idx0])) > 1 and np.amax(np.abs(lf[:,idx1])) > 1:
            SVD_grid[idx0,idx1] = S[1]/S[0]
        else:
            SVD_grid[idx0,idx1] = np.nan

plt.rcParams.update({'font.size': 12})

##### figure 4 in separability paper, SVD of 10 finger patches in M1 and S1
x0 = -.5
x1 = 9.5
y0 = -.5
y1 = 9.5

ix = np.tril_indices(10)
SVD_grid[ix] = np.nan

fig, ax = plt.subplots()
heatmap = ax.imshow(SVD_grid, origin='lower', extent=[x0,x1,y0,y1])
ax.set_title('SVD $\sigma_2 / \sigma_1$')
ax.set_xticks(np.arange(SVD_grid.shape[1]) , minor=False)
ax.set_yticks(np.arange(SVD_grid.shape[0]) , minor=False)
ax.invert_yaxis()
ax.xaxis.tick_top()
column_labels = ['','M1-d2','d3','d4','d5','S1-d1','d2','d3','d4','d5']
row_labels = ['M1-d1','M1-d2','M1-d3','M1-d4','M1-d5','S1-d1','S1-d2','S1-d3','S1-d4','S1-d5']
ax.set_xticklabels(column_labels, minor=False)
ax.set_yticklabels(row_labels, minor=False)
ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
fig.colorbar(heatmap)
heatmap.set_clim([0,0.8])

max_vals = np.amax(np.abs(lf),0)
max_vals = np.append(max_vals,[0,np.ceil(max_vals.max()/10)*10]) #so we get a gray scale between 0 and rounded max

colors = plt.cm.gray_r((max_vals-np.min(max_vals))/(np.max(max_vals)-np.min(max_vals)))

plt.figure()
img = plt.imshow(np.array([[max_vals.min(),max_vals.max()]]), cmap="Greys", origin='lower', extent=[x0,x1,y0,y1])
img.set_visible(False)
for i in range(10):
    if i != 0:
        plt.plot([i-.5,i+.5],[-.5,-.5], color = colors[i],linewidth=16)

    if i != 9:
        plt.plot([9.5,9.5],[i-.5,i+.5], color = colors[i],linewidth=16)
plt.colorbar(orientation="vertical")

plt.show()
