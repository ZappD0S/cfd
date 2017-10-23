import numpy as np
from itertools import product
import boundary_tools as bt
from utils import clockwise
import matplotlib.pyplot as plt
assert plt


def wrapping_rectangle(P, sizes):
    origin = [min(P[dim], P[dim] + sizes[dim]) for dim in range(2)]
    sizes = [abs(sizes[dim]) for dim in range(2)]
    return [(origin[0] - 0.5,            origin[1] - 0.5),
            (origin[0] + sizes[0] + 0.5, origin[1] - 0.5),
            (origin[0] + sizes[0] + 0.5, origin[1] + sizes[1] + 0.5),
            (origin[0] - 0.5,            origin[1] + sizes[1] + 0.5)]


class CellData:
    class WallSection():
        def __init__(self, midpt, normal, L):
            self.midpt = midpt
            self.normal = normal
            self.L = L

    no_free_perimeter = dict(((i, j), 0) for i, j in clockwise)
    default_perimeter = dict(((i, j), 2) for i, j in clockwise)

    def __init__(self):
        self.wall_sections = []
        self.free_perimeter = self.default_perimeter

    def add_wall_section(self, midpt, normal, L):
        self.wall_sections.append(self.WallSection(midpt, normal, L))


class GridData():
    default_cell_bounds = [[-1, 1], [-1, 1]]

    def __init__(self, shape, wall, dx, cell_bounds=None):

        self.ny, self.nx = self.shape = shape
        self.dx = dx

        if cell_bounds is None:
            self.cell_bounds = self.default_cell_bounds
        else:
            self.cell_bounds = cell_bounds

        self.walls = [bt.add_intersections(wall)]
        for y, sign in zip((0, self.ny - 1), (-1, 1)):
            side_wall = wrapping_rectangle((-1, y), (self.nx + 1, 3*sign))
            side_wall = bt.add_intersections(side_wall)
            self.walls.append(side_wall)

        self._centers_set = []
        for wall in self.walls:
            raw_centers = list(bt.data_cells(wall, self.cell_bounds))
            self._centers_set.append(
                [center for center in raw_centers if self.is_on_grid(center[::-1])])

        self._data_points_mask = self.build_mask_from_f(
            lambda ij: any(ij[::-1] in centers
                           for centers in self._centers_set))

        self._inside_wall_mask = np.any(
            [self.build_mask_from_f(
             lambda ij: bt.is_inside_wall(ij[::-1], wall))
             for wall in self.walls], axis=0)

        self._init_array()
        self.system_cells = self.build_mask_from_f(self.is_system)
        self.fluid_cells = self.build_mask_from_f(self.is_fluid)

        self._build_midpts_data()

    # Initialization functions

    def _init_array(self):
        self._array = np.empty(self.shape, dtype=object)
        self.max_n_ws = 0
        self.max_n_ws_midpts = 0

        for i in range(self.ny):
            for j in range(self.nx):
                self._array[i, j] = CellData()
                if self.is_inside_wall((i, j)) and not self.is_data_point((i, j)):
                    self._array[i, j].free_perimeter = CellData.no_free_perimeter

        for wall, centers in zip(self.walls, self._centers_set):
            for (j, i) in centers:
                subpts_set = bt.get_subpts((j, i), wall, self.cell_bounds)
                if subpts_set:
                    self._array[i, j].free_perimeter = \
                        bt.compute_free_perimeter(
                            subpts_set, (j, i), wall, self.cell_bounds)

                    n_subpts = len(subpts_set)
                    if n_subpts > self.max_n_ws:
                        self.max_n_ws = n_subpts

                    for subpts in subpts_set:
                        midpts, normals, L = bt.compute_len_midpt_normal(
                            subpts, (j, i), wall)
                        self._array[i, j].add_wall_section(
                            *midpts, *normals, L)

    def _build_midpts_data(self):
        self.midpts_set = [[] for _ in range(self.max_n_ws)]
        self.midpts_indices_set = []
        self.midpts_cells_set = [np.zeros(self.shape, dtype=np.bool)
                                 for _ in range(self.max_n_ws)]

        for i in range(self.ny):
            for j in range(self.nx):
                if self.is_data_point((i, j)) and self.is_system((i, j)):
                    for n, ws in enumerate(self.cell_ws_iterator((i, j))):
                        self.midpts_cells_set[n][i, j] = True
                        self.midpts_set[n].append(ws.midpt)

        for n in range(self.max_n_ws):
            self.midpts_set[n] = np.array(self.midpts_set[n])
            self.midpts_indices_set.append(
                self.build_indices(self.midpts_cells_set[n]))

    # Grid properties functions

    def is_on_grid(self, ij):
        return all(0 <= ij[dim] <= self.shape[dim]-1 for dim in range(2))

    def is_inside_wall(self, ij):
        if not self.is_on_grid(ij):
            ij = self.nearest_on_grid(ij)
        return self._inside_wall_mask[ij]

    def is_outside_wall(self, ij):
        return not self.is_inside_wall(ij)

    def is_data_point(self, ij):
        return self.is_on_grid(ij) and self._data_points_mask[ij]

    def is_outer_near_wall(self, ij):
        return self.is_data_point(ij) and self.is_outside_wall(ij)

    def is_fluid(self, ij):
        return self.is_outside_wall(ij) or self.is_data_point(ij)

    def is_system(self, ij):
        return (self.is_on_grid(ij) and
                self.is_fluid(ij) and
                not self.is_const_pressure(ij))

    def is_const_pressure(self, ij):
        return (0 <= ij[0] < self.shape[0] and
                (ij[1] <= 0 or ij[1] >= self.shape[1] - 1))

    def is_const_pressure_dir(self, ij1, ij2):
        return (self._is_near(self.is_const_pressure, ij1) and
                ij1[0] == ij2[0] and abs(ij2[1]-ij1[1]) == 1)

    # Grid data getter functions

    def cell_data_iterator(self, ij):
        if self.is_data_point(ij):
            for ws in self._array[ij].wall_sections:
                yield ws.midpt, ws.normal, ws.L
        else:
            raise ValueError('No data here..')

    def cell_ws_iterator(self, ij):
        if self.is_data_point(ij):
            for ws in self._array[ij].wall_sections:
                yield ws
        else:
            raise ValueError('No data here..')

    def free_perimeter(self, ij):
        return self._array[ij].free_perimeter.items()

    def n_indices(self, indices):
        return indices[indices >= 0][-1] + 1

    # Utilities
    # TODO: Should we move this in the utils.py file?

    def nearest_on_grid(self, ij):
        grid_ij = list(ij)
        for dim in range(2):
            if ij[dim] < 0:
                grid_ij[dim] = 0
            elif ij[dim] >= self.shape[dim]:
                grid_ij[dim] = self.shape[dim] - 1

        return tuple(grid_ij)

    def _is_near(self, f, ij):
        i, j = ij
        return not f(ij) and any(f((i+si, j+sj)) for si, sj in clockwise)

    def build_mask_from_f(self, f):
        mask = np.zeros(self.shape, dtype=np.bool)
        for i in range(self.ny):
            for j in range(self.nx):
                if f((i, j)):
                    mask[i, j] = True
        return mask

    def build_indices(self, *masks, op='or'):
        indices = np.full(self.shape, -1, dtype=np.int)
        count = 0
        for i in range(self.ny):
            for j in range(self.nx):
                if self._evaluate((i, j), masks, op):
                    indices[i, j] = count
                    count += 1

        return indices

    def _evaluate(self, ij, masks, op='or'):
        results = []
        for mask in masks:
            if hasattr(mask, '__call__'):
                results.append(mask(ij))
            elif hasattr(mask, '__getitem__'):
                if self.is_on_grid(ij):
                    results.append(mask[ij])
                else:
                    results.append(False)

        if op == 'or':
            return any(results)
        elif op == 'and':
            return all(results)
        else:
            raise ValueError('Invalid or unknown operation.')
