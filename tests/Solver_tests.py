import numpy as np
import matplotlib.pyplot as plt
import pickle
import context
from cfd.GridData import GridData
from cfd.Solver import Solver

assert context

# ny, nx = sizes = (80*2, 200*2)
ny, nx = sizes = (80, 200)

dx = 1.
dt = 0.02

vel0 = 140.

# x0, y0 = 0.2*sizes[1], sizes[0]/2
# R = 0.3*ny

x0, y0 = 0.3*sizes[1], sizes[0]/2
R = 0.18*sizes[0]

theta = np.linspace(0, 2*np.pi, 1000, endpoint=False)
x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)
circle = list(zip(x, y))
data = GridData(sizes, circle, dx)

# pickle.dump(data, open("data.p", "wb"))
# data = pickle.load(open("data.p", "rb"))

solv = Solver(data, dt, vel0=vel0, margin=0.1)

for i in range(400):
    solv.project()
    solv.advect()
    print(i)

plt.matshow(np.sqrt(solv._uv[0]**2 + solv._uv[1]**2))
plt.show()
