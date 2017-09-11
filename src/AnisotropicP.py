import numpy as np
from utils import ds_dict
import matplotlib.pyplot as plt
assert plt


class AnisotropicP:

    def __init__(self, p, uv, data):
        self.ny, self.nx = self.shape = data.shape
        self.dx = data.dx
        self.data = data

        self._uv = uv
        self._p = p

        self._inner_wall_mask = (data._inside_wall_mask & ~data._data_points_mask)
        self._build_offshore_data()

    def _build_offshore_data(self):
        self._offshore = self.data.system_cells.copy()
        self._offshore[0, :] = False
        self._offshore[-1, :] = False

        self._offshore_dir_masks = np.zeros((4, *self.shape), dtype=np.bool)
        self._bnd_dir_masks = np.zeros((4, *self.shape), dtype=np.bool)
        self._shifted_bnd_dir_masks = np.zeros((4, *self.shape), dtype=np.bool)
        offshore_opp_dir_mask = np.zeros(self.shape, dtype=np.bool)

        for (dim, sign), n in ds_dict.items():
            regular = self._shift_slice(dim, 0)
            forward = self._shift_slice(dim, sign)
            backward = self._shift_slice(dim, -sign)
            self._offshore_dir_masks[n][forward] = self._offshore[regular]

            offshore_opp_dir_mask[...] = 0
            offshore_opp_dir_mask[backward] = self._offshore[regular]
            self._bnd_dir_masks[n] = (offshore_opp_dir_mask & ~self._offshore)
            self._shifted_bnd_dir_masks[n][regular] = self._bnd_dir_masks[n][backward]

    def _shift_slice(self, dim, shift):
        #  in tutti gli array le dimensioni spaziali sono invertite.
        return tuple(slice(1+shift, self.shape[i]-1+shift)
                     if i != dim else slice(None) for i in range(2))

    def update_solver_bnd_p(self, bnd_p):
        bnd_p[...] = 0
        for (dim, sign), n in ds_dict.items():
            bnd_p[self._bnd_dir_masks[n]] += \
                sign*self.dx*self._uv[dim][self._bnd_dir_masks[n]]

    def update_vels(self):
        for dim in range(2):
            n_p = ds_dict[(dim, 1)]
            n_m = ds_dict[(dim, -1)]

            self._p[self._bnd_dir_masks[n_p]] = \
            self._p[self._shifted_bnd_dir_masks[n_p]]
            self._p[self._bnd_dir_masks[n_m]] = \
            self._p[self._shifted_bnd_dir_masks[n_m]]

            self._uv[dim][self._offshore] -= (
                self._p[self._offshore_dir_masks[n_p]]
                - self._p[self._offshore_dir_masks[n_m]])/(2*self.dx)

            self._uv[dim][self._inner_wall_mask] = 0
