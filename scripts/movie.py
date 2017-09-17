import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import context
from GridData import GridData
from Solver import Solver
# from matplotlib.animation import FuncAnimation, writers

assert context

# ffmpeg -framerate 60 -i frames/abs-%d.png -y abs.mp4
# ffmpeg -framerate 60 -i frames/curl-%d.png -y curl.mp4

ny, nx = sizes = (120, 360)
dpi = 100
dx = 1.
dt = 0.01

# x0, y0 = 0.2*sizes[1], sizes[0]/2
# R = 0.18*ny
#
# theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
# x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
# circle = list(zip(x, y))
#
# data = GridData(sizes, circle, dx)
# pickle.dump(data, open("data_large.p", "wb"))

data = pickle.load(open("data_large.p", "rb"))

N_per_vel = 1000
# max_vel = 200
# n_vels = 2
# vels = np.linspace(0, max_vel, n_vels+1)[1:]

vels = [50, 100, 200]
# vels = [100]
n_vels = len(vels)

solv = Solver(data, dt, vel0=0, margin=0.1, spline_type='bilinear')

figs = []
axs = []
ims = []
fig1, ax1 = plt.subplots()
ax1.axis('off')
figs.append(fig1)
axs.append(ax1)
im_abs = ax1.imshow(solv.get_magnitude_field(), interpolation='nearest')
ims.append(im_abs)

# fig2, ax2 = plt.subplots()
# ax2.axis('off')
# figs.append(fig2)
# axs.append(ax2)
# im_curl = ax2.imshow(solv.get_curl_field(), interpolation='nearest')
# ims.append(im_curl)

max_len = 0

for fig, ax, im in zip(figs, axs, ims):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad="1%")
    fig.colorbar(im, cax=cax)


for n, vel in enumerate(vels):
    solv.update_inlet_vels(vel)
    im_abs.set_norm(matplotlib.colors.Normalize(vmin=0, vmax= 2*vel))
    # im_curl.set_norm(matplotlib.colors.Normalize(vmin=-1.5*vel, vmax=1.5*vel))

    for i in range(N_per_vel):
        prn_str = "%.1f%%" % float(100 * (n*N_per_vel + i)/(n_vels*N_per_vel))
        if len(prn_str) > max_len:
            max_len = len(prn_str)
        print(prn_str + (max_len - len(prn_str))*' ', end='\r')
        solv.project()
        im_abs.set_data(solv.get_magnitude_field())
        # im_abs.autoscale()
        fig1.savefig(os.path.abspath('frames/abs-%d.png') % (n*N_per_vel + i),
                     bbox_inches='tight', dpi=120)

        # im_curl.set_data(solv.get_curl_field())
        # im_curl.autoscale()
        # fig2.savefig(os.path.abspath('frames/curl-%d.png') % (n*N_per_vel + i),
        #              bbox_inches='tight', dpi=120)
        solv.advect()

# ani = FuncAnimation(fig, animate, frames=N_per_vel*n_vels, interval=30)
# writer = writers['ffmpeg'](fps=30)
# ani.save('test.mp4', writer=writer, dpi=100,
#          savefig_kwargs={'bbox_inches': 'tight', 'pad_inches': 0})
