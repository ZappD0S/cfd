import numpy as np
# from scipy.signal import convolve2d

neighbor_mean = np.array([[2**-0.5, 1, 2**-0.5],
                          [1,       1,       1],
                          [2**-0.5, 1, 2**-0.5]])
neighbor_mean /= np.sum(neighbor_mean)


def advect(u, v, dx, dt):

    u0, v0 = u.copy(), v.copy()
    u_back, v_back = np.zeros_like(u), np.zeros_like(u)
    u_forw, v_forw = np.zeros_like(u), np.zeros_like(u)

    def material_back2(x, x0, u, v, k):
        x[1:-1, 1:] -= k*u[1:-1, 1:]*(x0[1:-1, 1:] - x0[1:-1, :-1])
        x[1:, 1:-1] -= k*v[1:, 1:-1]*(x0[1:, 1:-1] - x0[:-1, 1:-1])

    def material_forw2(x, x0, u, v, k):
        x[1:-1, :-1] -= k*u[1:-1, :-1]*(x0[1:-1, 1:] - x0[1:-1, :-1])
        x[:-1, 1:-1] -= k*v[:-1, 1:-1]*(x0[1:, 1:-1] - x0[:-1, 1:-1])
    #
    # def material_back(x, u, v):
    #     return -((u[1:-1, 1:-1] + u[1:-1, :-2])/2
    #              * (x[1:-1, 1:-1] - x[1:-1, :-2])/dx
    #              + (v[1:-1, 1:-1] + v[:-2, 1:-1])/2
    #              * (x[1:-1, 1:-1] - x[:-2, 1:-1])/dx)
    #
    # def material_forw(x, u, v):
    #     return -((u[1:-1, 2:] + u[1:-1, 1:-1])/2
    #              * (x[1:-1, 2:] - x[1:-1, 1:-1])/dx
    #              + (v[2:, 1:-1] + v[1:-1, 1:-1])/2
    #              * (x[2:, 1:-1] - x[1:-1, 1:-1])/dx)

    # def material_back(x, u, v):
    #     return -((u[1:-1, 1:-1]*x[1:-1, 1:-1] - u[1:-1, :-2]*x[1:-1, :-2])/dx
    #              + (v[1:-1, 1:-1]*x[1:-1, 1:-1] - v[:-2, 1:-1]*x[:-2, 1:-1])/dx)
    #
    # def material_forw(x, u, v):
    #     return -((u[1:-1, 2:]*x[1:-1, 2:] - u[1:-1, 1:-1]*x[1:-1, 1:-1])/dx
    #              + (v[2:, 1:-1]*x[2:, 1:-1] - v[1:-1, 1:-1]*x[1:-1, 1:-1])/dx)

    def material_back(x, x0, u, v, k):
        x[1:-1, 1:] -= k*(u[1:-1, 1:]*x0[1:-1, 1:]
                          - u[1:-1, :-1]*x0[1:-1, :-1])
        x[1:, 1:-1] -= k*(v[1:, 1:-1]*x0[1:, 1:-1]
                          - v[:-1, 1:-1]*x0[:-1, 1:-1])

    def material_forw(x, x0, u, v, k):
        x[1:-1, :-1] -= k*(u[1:-1, 1:]*x0[1:-1, 1:]
                           - u[1:-1, :-1]*x0[1:-1, :-1])
        x[:-1, 1:-1] -= k*(v[1:, 1:-1]*x0[1:, 1:-1]
                           - v[:-1, 1:-1]*x0[:-1, 1:-1])

    k = dt/dx

    u_back[...] = u0
    v_back[...] = v0
    material_forw2(u_back, u0, u0, v0, k)
    material_forw2(v_back, v0, u0, v0, k)

    # questi rendono le derivate ai margini nulle per material_forw
    # u_forw[1:-1, -1] = u_forw[1:-1, -2]
    # u_forw[-1, 1:-1] = u_forw[-2, 1:-1]
    #
    # v_forw[1:-1, -1] = v_forw[1:-1, -2]
    # v_forw[-1, 1:-1] = v_forw[-2, 1:-1]

    # u_forw[1:-1, 0] = u_forw[1:-1, 1]
    # u_forw[0, 1:-1] = u_forw[1, 1:-1]
    #
    # v_forw[1:-1, 0] = v_forw[1:-1, 1]
    # v_forw[0, 1:-1] = v_forw[1, 1:-1]

    # u_back[1:-1, 1:-1] = \
    #     u_forw[1:-1, 1:-1] - dt*material_forw(u_forw, u_forw, v_forw)
    # v_back[1:-1, 1:-1] = \
    #     v_forw[1:-1, 1:-1] - dt*material_forw(v_forw, u_forw, v_forw)
    #
    # u_corr = u0 - (u_back-u0)/2
    # v_corr = v0 - (v_back-v0)/2

    # u[1:-1, 1:-1] = u_corr[1:-1, 1:-1] + dt*material_back(u_corr, u_corr, v_corr)
    # v[1:-1, 1:-1] = v_corr[1:-1, 1:-1] + dt*material_back(v_corr, u_corr, v_corr)

    u[...] = (u0 + u_back)/2
    v[...] = (v0 + v_back)/2

    material_back2(u, u_back, u_back, v_back, k/2)
    material_back2(v, v_back, u_back, v_back, k/2)
