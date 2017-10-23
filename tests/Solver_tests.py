import numpy as np
import matplotlib.pyplot as plt
import pickle
import context
from GridData import GridData
from Solver import Solver


#ny, nx = sizes = (80*4, 200*4)
ny, nx = sizes = (80, 200)
#ny, nx = sizes = (120, 360)

dx = 1.

############################################################

#x0, y0 = 0.3*sizes[1], sizes[0]/2
#R = 0.2*ny

x0, y0 = 0.2*sizes[1], sizes[0]/2
R = 0.18*ny


theta = np.linspace(np.pi/2, 3/2 * np.pi, 1000, endpoint=False)
x, y = x0+R*np.cos(theta), y0+R*np.sin(theta)

P0 = [x[-1], y[-1]]

x = np.append(x, P0[0] + 4*R)
y = np.append(y, P0[1])
figure = list(zip(x, y))

# circle = list(zip(x, y))
data = GridData(sizes, figure, dx)
# pickle.dump(data, open("./data.p", "wb"))

############################################################

dt = 0.01
vel0 = 150.

# data = pickle.load(open("./data.p", "rb"))
solv = Solver(data, dt, vel0=vel0, margin=0.1, spline_type="bilinear")

for i in range(500):
    solv.project()
    solv.advect()
    print(i)

solv.close()
plt.matshow(solv.get_magnitude_field())
plt.title("vel0: %.1e, dt: %.e" % (vel0, dt))
plt.show()
