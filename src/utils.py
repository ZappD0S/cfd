import numpy as np
from itertools import tee, chain
from ctypes import POINTER, pointer, c_double

clockwise = [(0, 1), (1, 0), (0, -1), (-1, 0)]

ds_dict = {}
for n in range(4):
    sign = sum(clockwise[n])
    dim = clockwise[n].index(0)
    ds_dict[(dim, sign)] = n


def pairs(iterable, last=False):
    if last:
        it = iter(iterable)
        first = [next(it)]
        iterable = chain(first, it, first)

    its = tee(iterable, 2)
    for i, it in enumerate(its):
        for _ in range(i):
            next(it, None)

    return zip(*its)


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

# Ctypes stuff

def build_c_vecs(vecs):
    length = vecs[0].size if type(vecs[0]) == np.ndarray else len(vecs[0])
    vec_type = c_double * length
    c_vecs = (POINTER(vec_type) * len(vecs))()
    for n, vec in enumerate(vecs):
        if type(vec) == np.ndarray:
            c_vecs[n] = vec.ctypes.data_as(POINTER(vec_type))
        else:
            c_vecs[n] = pointer(vec_type(*vec))

    return c_vecs
