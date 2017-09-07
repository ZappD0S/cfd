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
        self._bnd_p = np.zeros((2, *p.shape), dtype=np.double)

        self._bnd_dir_masks = data.build_dir_masks(data.is_fluid)

        self._inside_wall_mask = data.build_mask_from_f(data.is_inside_wall)

        # TODO: le multiplicities devono essere contante solo per dimensioni
        # uguali! Comunque adesso non servono
        # self._multiplicities = np.sum(self._bnd_masks, axis=-1)

        self._build_offshore_data()

        self._nw_dir_masks = data.build_dir_masks(
            self.data.is_outside_wall)
        self._shifted_nw_dir_masks = self._build_shifted_masks(
            self._nw_dir_masks)

    def _build_offshore_data(self):
        bnd_cells = np.any(self._bnd_dir_masks, axis=0)
        self._offshore = self.data.system_cells & ~bnd_cells
        self._offshore_dir_masks = np.zeros((4, *self.shape), dtype=np.bool)

        for (dim, sign), n in ds_dict.items():
            regular = self._shift_slice(dim, 0)
            shifted = self._shift_slice(dim, sign)
            self._offshore_dir_masks[n][shifted] = self._offshore[regular]

    def _build_shifted_masks(self, masks):
        shifted_masks = np.zeros_like(masks)
        for dim in range(2):
            regular = self._shift_slice(dim, 0)
            for sign in (1, -1):
                n = ds_dict[(dim, sign)]
                shifted = self._shift_slice(dim, sign)
                shifted_masks[n][shifted] = masks[n][regular]
        return shifted_masks

    def _shift_slice(self, dim, shift):
        #  in tutti gli array le dimensioni spaziali sono invertite.
        return tuple(slice(1+shift, self.shape[i]-1+shift)
                     if i != dim else slice(None) for i in range(2))

    def compute_solver_bnd_p(self, bnd_p):
        bnd_p[...] = 0
        for (dim, sign), n in ds_dict.items():
            bnd_p[self._bnd_dir_masks[n]] += \
                sign*self.dx*self._uv[dim][self._bnd_dir_masks[n]]

    def update_vels(self):
        for dim in range(2):
            n_p = ds_dict[(dim, 1)]
            n_m = ds_dict[(dim, -1)]

            self._uv[dim][self._offshore] -= (
                self._p[self._offshore_dir_masks[n_p]]
                - self._p[self._offshore_dir_masks[n_m]])/(2*self.dx)

            self._uv[dim][self._inside_wall_mask] = 0
