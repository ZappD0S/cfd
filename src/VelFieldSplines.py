import numpy as np
from ctypes import CDLL, POINTER, c_double, c_ulong, c_int, Structure
from utils import build_c_vecs


class Spline2DGroup(Structure):
    pass


splines_lib = CDLL('../bin/splines.so')
splines_lib.Spline2DGroup_alloc.restype = POINTER(Spline2DGroup)

advection_lib = CDLL('../bin/advection.so')


class VelFieldSplines:
    def __init__(self, shape, dt):
        self.shape = (c_ulong*2)(*shape)
        self.dt = c_double(dt)
        self.spl_g = splines_lib.Spline2DGroup_alloc(self.shape, c_int(2))

    def init_data(self, wallvels_set, midpts_set):
        self.midpts_set = []
        self.wallvels_set = []

        for midpts, wallvels in zip(midpts_set, wallvels_set):
            self.midpts_set.append(build_c_vecs(midpts))
            self.wallvels_set.append(build_c_vecs(wallvels))
        return self

    def init_vels(self, uv):
        self.c_uv = build_c_vecs(uv)
        self.update_splines()
        return self

    def update_splines(self):
        splines_lib.Spline2DGroup_set(self.spl_g, self.c_uv)

    def update_wallvels(self):
        self.update_splines()
        for midpts, wallvels in zip(self.midpts_set, self.wallvels_set):
            splines_lib.Spline2DGroup_array_eval(
                self.spl_g, midpts, wallvels, c_ulong(len(midpts)))

    def advect_vels(self):
        self.update_splines()
        advection_lib.advect_vels(self.spl_g, self.c_uv, self.dt)

    def close(self):
        splines_lib.Spline2DGroup_free(self.spl_g)
