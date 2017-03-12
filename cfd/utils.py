import numpy as np

clockwise = [(0, 1), (1, 0), (0, -1), (-1, 0)]

ds_dict = {}
for n in range(4):
    sign = sum(clockwise[n])
    dim = clockwise[n].index(0)
    ds_dict[(dim, sign)] = n


def spline_coeffs(center, mid_pt):
    ds = [b-a for a, b in zip(center, reversed(mid_pt))]
    return np.outer(*[[1-ds[i], ds[i]] for i in range(2)])


def wrapping_rectangle(P, sizes):
    origin = [min(P[dim], P[dim] + sizes[dim]) for dim in range(2)]
    sizes = [abs(sizes[dim]) for dim in range(2)]
    return [(origin[0] - 0.5,            origin[1] - 0.5),
            (origin[0] + sizes[0] + 0.5, origin[1] - 0.5),
            (origin[0] + sizes[0] + 0.5, origin[1] + sizes[1] + 0.5),
            (origin[0] - 0.5,            origin[1] + sizes[1] + 0.5)]


def build_inlet_vels(shape, margin):
    def smooth_step(x):
        return 6*x**5 - 15*x**4 + 10*x**3

    ny, nx = shape
    vels = np.ones(ny)
    n_margin = round(ny*margin)
    edge_vels = smooth_step(np.linspace(0, 1, n_margin))
    vels[:n_margin] = edge_vels
    vels[-n_margin:] = edge_vels[::-1]
    return vels


def get_time_step(u, v, dx):
    vel_max = ((u**2 + v**2)**0.5).max()
    print('vel_max:', vel_max)
    if vel_max:
        return 1 * dx/vel_max
    else:
        raise ValueError('u and v arrays have all zero values.')
