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

# ffmpeg -r 60 -i frames/abs-%d.png -y abs.mp4
# ffmpeg -r 60 -i frames/curl-%d.png -y curl.mp4

ny, nx = sizes = (100, 300)
dpi = 100
dx = 1.
dt = 0.01

# x0, y0 = 0.2*sizes[1], sizes[0]/2
# R = 0.22*ny
#
# theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
# x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
# circle = list(zip(x, y))
#
# data = GridData(sizes, circle, dx)
# pickle.dump(data, open("data.p", "wb"))

data = pickle.load(open("data.p", "rb"))

max_vel = 200
N_per_vel = 1000
# n_vels = 2
# vels = np.linspace(0, max_vel, n_vels+1)[1:]

vels = [50, 100, 200]
n_vels = len(vels)


solv = Solver(data, dt, vel0=0, margin=0.1)
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
ax1.axis('off')
ax2.axis('off')

im_abs = ax1.imshow(solv.get_magnitude_field(), interpolation='nearest')
im_curl = ax2.imshow(solv.get_curl_field(), interpolation='nearest')


for fig, ax, im in zip((fig1, fig2), (ax1, ax2), (im_abs, im_curl)):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad="1%")
    fig.colorbar(im, cax=cax)


for n, vel in enumerate(vels):
    solv.update_inlet_vels(vel)
    im_abs.set_norm(matplotlib.colors.Normalize(vmin=0, vmax=3*vel))
    im_curl.set_norm(matplotlib.colors.Normalize(vmin=-1.5*vel, vmax=1.5*vel))

    for i in range(N_per_vel):
        print("%.1f%%" % float(100 * (n*N_per_vel + i)/(n_vels*N_per_vel)))
        solv.project()
        im_abs.set_data(solv.get_magnitude_field())
        # im_abs.autoscale()
        im_curl.set_data(solv.get_curl_field())
        # im_curl.autoscale()

        fig1.savefig(os.path.abspath('frames/abs-%d.png') % (n*N_per_vel + i),
                     bbox_inches='tight', dpi=120)
        fig2.savefig(os.path.abspath('frames/curl-%d.png') % (n*N_per_vel + i),
                     bbox_inches='tight', dpi=120)

        solv.advect()


# print(sum(sizes[0] for sizes in divider.get_horizontal_sizes(im)))
# fig.subplots_adjust(left=0.05, bottom=0, right=0.9, top=1, wspace=0, hspace=0)
#
# plt.show()
# sizes = fig.get_size_inches()
# print(sizes)
#
# extent = im.get_extent()
# print(extent)
# fig.set_size_inches((sizes[0], sizes[0]*ny/nx))
# sizes = fig.get_size_inches()
# print(sizes)
# plt.show()


# last_change = None
#
#
# def animate(n):
#     print('n:', n)
#     global last_change
#     if n % N_per_vel == 0 and n != last_change:
#         last_change = n
#         vel = vels[n//N_per_vel % n_vels]
#         print('vel', vel)
#         solv.update_inlet_vels(vel)
#     solv.project()
#     im.set_data(solv.get_magnitude_field())
#     im.autoscale()
#     solv.advect()
#     return [im]

# ani = FuncAnimation(fig, animate, frames=N_per_vel*n_vels, interval=30)
# writer = writers['ffmpeg'](fps=30)
# ani.save('test.mp4', writer=writer, dpi=100,
#          savefig_kwargs={'bbox_inches': 'tight', 'pad_inches': 0})
