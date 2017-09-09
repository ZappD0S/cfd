import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized
from utils import build_inlet_vels, clockwise, ds_dict
# from VelFieldSplines_old import VelFieldSplines
from VelFieldSplines import VelFieldSplines

from AnisotropicP import AnisotropicP
import matplotlib.pyplot as plt
assert plt

# TODO: rename FluidFlow?
class Solver:

    def __init__(self, data, dt, vel0, margin=0.2,
                 decay_magnitude=1e-6, spline_type = "bilinear"):
        self.ny, self.nx = self.shape = data.shape
        self.data = data
        self.dx = data.dx
        self.dt = dt
        self.inlet_vel = vel0
        self.margin = margin
        self.decay_magnitude = decay_magnitude

        self._uv = np.zeros((2, *self.shape), dtype=np.double)
        self._default_inlet_vels = build_inlet_vels(self.shape, margin)
        self._init_inlet_vels()

        self._wallvels_set = [np.empty((2, len(midpts)), dtype=np.double)
                              for midpts in self.data.midpts_set]

        self.vfsplines = VelFieldSplines(
            self.shape, dt/self.dx, spline_type=spline_type)
        self.vfsplines.init_vels(self._uv).init_data(
            self._wallvels_set, self.data.midpts_set)

        self._system_indices = data.build_indices(data.is_system)
        self._n_sc = data.n_indices(self._system_indices)

        self._fluid_indices = data.build_indices(data.is_fluid)
        self._n_fc = data.n_indices(self._fluid_indices)

        self._coeffs = np.zeros(self._n_sc)
        self.div = np.zeros(self.shape)

        self._build_div_matrices()
        self._build_wall_div_matrices()
        self._update_div()

        self._p = np.zeros(self.shape, dtype=np.double)
        self._bnd_p = np.zeros_like(self._p)
        self.anis_p = AnisotropicP(self._p, self._uv, data)

        self._build_p_matrix()
        self._LUsolve = factorized(self._p_mat)

    # TODO: in the building functions, should we use dicts to store data?
    # Initialization functions

    def _build_p_matrix(self):
        rows, cols, vals = [], [], []

        for i in range(self.ny):
            for j in range(self.nx):
                if self.data.is_system((i, j)):
                    row = self._system_indices[i, j]
                    val = -4

                    for (si, sj), (dim, sign) in zip(clockwise, ds_dict):
                        if self.data.is_system((i+si, j+sj)):
                            rows.append(row)
                            cols.append(self._system_indices[i+si, j+sj])
                            vals.append(1)
                        else:
                            val += 1

                    rows.append(row)
                    cols.append(row)
                    vals.append(val)

        self._p_mat = coo_matrix((vals, (rows, cols)),
                                 shape=(self._n_sc, self._n_sc)).tocsc()

    def _build_div_matrices(self):
        row_set = [[], []]
        col_set = [[], []]
        val_set = [[], []]

        for i in range(self.ny):
            for j in range(self.nx):
                if self.data.is_system((i, j)):
                    row = self._system_indices[i, j]
                    for (si, sj), fp in self.data.free_perimeter((i, j)):
                        if (self.data.is_fluid((i+si, j+sj)) and
                           not self.data.is_const_pressure_dir(
                                (i, j), (i+si, j+sj))):
                            dim = 1 if si else 0
                            sign = si + sj
                            row_set[dim].append(row)
                            col_set[dim].append(
                                self._fluid_indices[i+si, j+sj])
                            val_set[dim].append(sign*fp/4)

        self._div_mats = [coo_matrix(
                            (val_set[dim], (row_set[dim], col_set[dim])),
                            shape=(self._n_sc, self._n_fc)).tocsc()
                          for dim in range(2)]

    def _build_wall_div_matrices(self):
        def empty():
            return [[] for _ in range(self.data.max_n_ws)]

        self._wall_div_mats = empty()

        for dim in range(2):
            row_set = empty()
            col_set = empty()
            val_set = empty()

            for i in range(self.ny):
                for j in range(self.nx):
                    if (self.data.is_data_point((i, j)) and
                       self.data.is_system((i, j))):
                        for n, ws in enumerate(
                                self.data.cell_ws_iterator((i, j))):
                            index = self.data.midpts_indices_set[n][i, j]
                            row_set[n].append(index)
                            col_set[n].append(index)
                            val_set[n].append(ws.L*ws.normal[dim]/4)

            for n, ws_mats in enumerate(self._wall_div_mats):
                size = len(val_set[n])
                ws_mats.append(coo_matrix(
                    (val_set[n], (row_set[n], col_set[n])),
                    shape=(size, size)).tocsc())

    def _init_inlet_vels(self):
        tau = -self.nx/(20*np.log(self.decay_magnitude))
        n = 0
        vel = self.inlet_vel
        while True:
            self._uv[0][:, n] = vel*self._default_inlet_vels
            n += 1
            vel = self.inlet_vel*np.exp(-n/tau)
            if vel < self.decay_magnitude:
                break

    # Update functions

    def update_inlet_vels(self, vel):
        if vel and not self.inlet_vel:
            self._init_inlet_vels()
        else:
            self._uv[0][:, 0] = vel*self._default_inlet_vels
        self.inlet_vel = vel

    def _update_div(self):
        self.div[self.data.system_cells] = sum(self._div_mats[dim].dot(
                                self._uv[dim][self.data.fluid_cells])
                            for dim in range(2))

        self.vfsplines.update_wallvels()
        for dim in range(2):
            for n in range(self.data.max_n_ws):
                self.div[self.data.midpts_cells_set[n]] += \
                    self._wall_div_mats[n][dim].dot(
                        self._wallvels_set[n][dim])

    def project(self):
        self.anis_p.compute_solver_bnd_p(self._bnd_p)
        self._coeffs[...] = (self.div[self.data.system_cells]*self.dx
                             - self._bnd_p[self.data.system_cells])

        self._p[self.data.system_cells] = self._LUsolve(self._coeffs)
        self._p[:, 0] = self._p[:, 1]
        self._p[:, -1] = self._p[:, -2]

        self.anis_p.update_vels()
        self._update_div()

    def advect(self):
        self.vfsplines.advect_vels()
        self._update_div()

    def close(self):
        self.vfsplines.close()

    # Getter functions

    def get_magnitude_field(self):
        return np.sqrt(self._uv[0]**2 + self._uv[1]**2)[1:-1, 1:-1]

    def get_curl_field(self):
        return (np.gradient(self._uv[1], self.dx, axis=1)
                - np.gradient(self._uv[0], self.dx, axis=0))[1:-1, 1:-1]
