import numpy as np
import matplotlib.pyplot as plt
import context
from GridData import GridData
assert context


def plot_rect(i, j, dir_):
    plt.plot([j, j, j+dir_[1], j+dir_[1], j],
             [i, i+dir_[0], i+dir_[0], i, i], 'green')

dx = 1

# ny, nx = sizes = (80, 200)


# x0, y0 = 0.3*sizes[1], sizes[0]/2
# theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
# R = 0.18*sizes[0]
# x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)

# ny, nx = sizes = (80, 200)
# x0, y0 = 0.3*sizes[1], sizes[0]/2
# R = 0.2*ny

ny, nx = sizes = (80, 200)
x0, y0 = 0.2*sizes[1], sizes[0]/2
R = 0.18*ny

theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
circle = list(zip(x, y))
data = GridData(sizes, circle, dx)

# plt.matshow(data._inside_wall_mask)
# plt.show()


for wall in data.walls:
    xs = [P[0] for P in wall]
    ys = [P[1] for P in wall]
    plt.plot(xs, ys, 'lightblue')

for i in range(data.ny):
    for j in range(data.nx):
        if data.is_data_point((i, j)):
            plt.plot(j, i, 'b.')
            for mid_pt, normal, L in data.cell_data_iterator((i, j)):
                xs, ys = [[mid_pt[i], mid_pt[i]+normal[i]] for i in range(2)]
                plt.plot(xs, ys, 'purple')
                plt.plot(*mid_pt, 'ro')
        if data._inside_wall_mask[i, j]:
            plt.plot(j, i, 'go')


plt.axis('equal')
plt.show()
