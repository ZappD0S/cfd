import numpy as np
from cfd.C.boundary_tools_utils import build_c_pts
from ctypes import c_double, c_ulong, POINTER, pointer, Structure, CDLL
import os


pwd = os.path.dirname(__file__)
spline_group_lib = CDLL(os.path.join(pwd, 'spline_group.so'))
advection_lib = CDLL(os.path.join(pwd, 'advection.so'))


class spline_group(Structure):
    pass


spline_group_lib.alloc_spline_group.restype = POINTER(spline_group)
spline_group_lib.spline_group_eval.restype = c_double


def build_c_uv(uv):
    ny, nx = uv[0].shape
    c_uv = (POINTER((c_double*nx)*ny)*2)()
    for dim in range(2):
        c_uv[dim] = pointer(np.ctypeslib.as_ctypes(uv[dim]))
    return c_uv


class VelFieldSpline:
    def __init__(self, shape, dt):
        ny, nx = self.shape = shape
        self.dt = c_double(dt)
        c_ny = c_ulong(ny)
        c_nx = c_ulong(nx)
        self.spl_g_set = (POINTER(spline_group)*2)()
        for dim in range(2):
            self.spl_g_set[dim] = spline_group_lib.alloc_spline_group(
                c_ny, c_nx)

    def init_wall_vels(self, wall_vels_set, midpts_set):
        self.c_midpts_set = []
        self.c_wv_set = []

        for uv_wvs, midpts in zip(wall_vels_set, midpts_set):
            self.c_midpts_set.append(build_c_pts(midpts))
            self.c_wv_set.append([np.ctypeslib.as_ctypes(uv_wvs[dim])
                                  for dim in range(2)])
        return self

    def init_vels(self, uv):
        self.c_uv = build_c_uv(uv)
        self._update_splines()
        return self

    def _update_splines(self):
        for dim in range(2):
            spline_group_lib.init_spline_group(
                self.spl_g_set[dim], self.c_uv[dim])

    def close(self):
        for spl_g in self.spl_g_set:
            spline_group_lib.free_spline_group(spl_g)
        self.__delattr__('spl_g_set')

    def update_wall_vels(self):
        self._update_splines()
        for uv_wvs, midpts in zip(self.c_wv_set, self.c_midpts_set):
            for dim, wvs in enumerate(uv_wvs):
                spline_group_lib.spline_group_array_eval(
                    self.spl_g_set[dim], wvs, midpts, c_ulong(len(wvs)))

    def advect_vels(self):
        self._update_splines()
        advection_lib.advect_vels(self.spl_g_set, self.c_uv, self.dt)
