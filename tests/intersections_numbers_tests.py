import numpy as np
import context
import cfd.boundary_tools as bt
from cfd.C.boundary_tools_utils import intersections_number
import matplotlib.pyplot as plt

assert context

ny, nx = sizes = (100, 300)
x0, y0 = 0.2*sizes[1], sizes[0]/2
theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
R = 0.22*ny
x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
circle = list(zip(x, y))

circle = bt.add_intersections(circle)

normal = [0.68225361, 0.73111559]
P0 = (44.970860153335522, 33.93385858964718)
print(intersections_number(normal, P0, circle))

fig, ax = plt.subplots()

xs = [P[0] for P in circle]
ys = [P[1] for P in circle]
ax.plot(xs, ys)

PS = [(30.029140, 32.159475), (29.988114, 32.197759),
      (-0.030896, 0.028832), (0.029140, -0.027192)]

for P in PS:
    ax.plot(*[P0[i]+P[i] for i in range(2)], 'go')

ax.plot(*P0, 'ro')
xs, ys = [[P0[i], P0[i]+normal[i]] for i in range(2)]
ax.plot(xs, ys, 'purple')
plt.axis('equal')
plt.show()
