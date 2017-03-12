import numpy as np
import matplotlib.pyplot as plt
import context
from cfd.boundary_tools import get_subpts, occupied_cells, \
                               insert_intersections, compute_area, \
                               compute_free_perimeter, \
                               compute_len_midpt_normal
from libs.pairs import pairs

assert context
assert plt

ny, nx = 100, 200
sizes = (ny, nx)
dx = 0.1

x0, y0 = 50, 50
theta = np.linspace(0, 2*np.pi, 200, endpoint=False)
R = 35


# ny, nx = 100, 105
# sizes = (ny, nx)
# dx = 1
#
# x0, y0 = 50, 50
# theta = np.linspace(0, 2*np.pi, 100, endpoint=False)
# R = 15
x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
circle = list(zip(x, y))

circle = insert_intersections(circle)
centers = list(occupied_cells(circle))


dirs = [(1, 1), (-1, 1), (-1, -1), (1, -1)]

for i in [84]:
    # print("\n")
    for j in [42]:
        if (j, i) in centers:
            subpts_set = get_subpts((j, i), circle)
            vertexes = [[[c+s for c, s in zip((j, i), d)] for d in pair]
                        for pair in pairs(dirs, last=True)]
            A = compute_area(subpts_set, (j, i), circle)
            print('center:', (j, i), 'area:', A)
            fp = compute_free_perimeter(subpts_set, (j, i), circle)
            print(fp)

            plt.gca().set_aspect('equal', adjustable='box')
            for pair in vertexes:
                xs = [P[0] for P in pair]
                ys = [P[1] for P in pair]
                plt.plot(xs, ys, 'lightblue')

            for P in circle:
                plt.plot(*P, 'bo')
            for subpts in subpts_set:
                for P in subpts:
                    plt.plot(*P, 'ro')

            for subpts in subpts_set:
                L, mid_pt, normal = compute_len_midpt_normal(subpts, (j, i),
                                                             circle)
                print('L:', L)
                xs, ys = [[mid_pt[i], mid_pt[i]+normal[i]] for i in range(2)]
                plt.plot(xs, ys, 'purple')
            plt.show()
