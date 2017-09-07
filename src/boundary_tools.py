import numpy as np
from math import ceil, floor
from itertools import product
from functools import cmp_to_key
from ctypes import CDLL, c_double, c_int, c_ulong
from utils import pairs, build_c_vecs

inters_num_lib = CDLL('../bin/intersections_number.so')

clockwise_inds = [(0, 0), (0, 1), (1, 1), (1, 0)]
clockwise_ds = [(0, -1), (1, 1), (0, 1), (1, -1)]
signs = {(0, -1): 1, (1, 1): 1, (0, 1): -1, (1, -1): -1}


def add_intersections(pts):
    def find_y_with_x(P1, P2, x0):
        return (P2[1]-P1[1])/(P2[0]-P1[0]) * (x0-P1[0])+P1[1]

    def find_x_with_y(P1, P2, y0):
        return (P2[0]-P1[0])/(P2[1]-P1[1]) * (y0-P1[1])+P1[0]

    find_other_coord = (find_y_with_x, find_x_with_y)

    pts_w_inters = []
    added = 0
    checked = pts.copy()
    for index, (P1, P2) in enumerate(pairs(pts, last=True)):
        inters = []
        for dim in range(2):
            upper = max(P1[dim], P2[dim])
            lower = min(P1[dim], P2[dim])
            delta = ceil(upper)-floor(lower)

            if delta >= 2:
                # TODO: non mi piace j, troviamo un altro nome..
                for j in range(floor(lower)+1, ceil(upper)):
                    inter = np.empty(2)
                    inter[dim] = j
                    inter[1-dim] = find_other_coord[dim](P1, P2, j)
                    inter = tuple(inter)
                    if not any(np.allclose(inter, P) for P in checked):
                        inters.append(inter)
                        checked.append(inter)
                    # if inter not in inters:
                    #     inters.append(tuple(inter))

        inters = sort(P1, inters)
        pts_w_inters = (pts_w_inters[:added+index+1]
                        + inters+pts_w_inters[added+index+1:])
        added += len(inters)

    return pts_w_inters


def distance(P1, P2):
    return ((P2[0]-P1[0])**2+(P2[1]-P1[1])**2)**(1/2)


def sort(P0, pts):
    sorted_pts = [P0]
    while pts:
        remaining = [pts[i] for i in range(len(pts))]
        nearest = min(remaining, key=lambda P: distance(sorted_pts[-1], P))
        sorted_pts.append(nearest)
        pts.remove(nearest)
    return sorted_pts


def compute_len_midpt_normal(subpts, center, pts, N=1):
    def midpoint(P1, P2, t):
        return tuple(P1[dim]*(1-t)+P2[dim]*t for dim in range(2))

    def slope_ang(P1, P2):
        return np.arctan2(P2[1]-P1[1], P2[0]-P1[0])

    # quando i==1, subpts[:1] contiene solo un elemento e
    # pairs restituisce None, quindi sum([]) == 0
    partial_lengths = [sum(distance(p1, p2) for p1, p2 in pairs(subpts[:i]))
                       for i in range(1, len(subpts)+1)]
    L = partial_lengths[-1]
    mid_pts = []
    normals = []

    for n in range(N):
        L_frac = L * (n+1)/(N+1)
        if any(np.isclose(L_frac, l) for l in partial_lengths):
            nb_pts = next([subpts[i+s] for s in (-1, 0, 1)]
                          for i in range(1, len(subpts)-1)
                          if np.isclose(L_frac, partial_lengths[i]))
            mid_pts.append(nb_pts[1])
            theta = sum(slope_ang(p1, p2) for p1, p2 in pairs(nb_pts))/2
        else:
            nb_lengths = next((l1, l2) for l1, l2 in pairs(partial_lengths)
                              if l1 < L_frac < l2)
            nb_pts = next(subpts[i:i+2]
                          for i, (l1, l2) in enumerate(pairs(partial_lengths))
                          if l1 < L_frac < l2)
            mid_pt = midpoint(*nb_pts,
                              (L_frac-nb_lengths[0])
                              / (nb_lengths[1]-nb_lengths[0]))
            mid_pts.append(mid_pt)
            theta = slope_ang(*nb_pts)

        normal = np.array([np.cos(theta+np.pi/2), np.sin(theta+np.pi/2)])
        normal[np.abs(normal) < 1e-15] = 0
        normal *= get_inward_normal_sign(normal, mid_pts[-1], pts)
        normals.append(normal)

    return mid_pts, normals, L


def compute_free_perimeter(subpts_set, center, pts, cell_bounds):

    def get_vertices(cell_bounds, inds=[], n=0):
        if n < len(cell_bounds):
            return [get_vertices(cell_bounds, inds+[i], n+1)
                    for i in cell_bounds[n]]
        else:
            return tuple(inds)

    verts = get_vertices(cell_bounds)
    clockwise_verts = [verts[i][j] for i, j in clockwise_inds]
    cell_vertices = [tuple(c+s for c, s in zip(center, vertex))
                     for vertex in clockwise_verts]

    free_perimeter = dict(((i, j), 0) for i in (1, 0, -1) for j in (1, 0, -1)
                          if abs(i) != abs(j))

    def get_common_dim_sign(P1, P2):
        for dim in range(2):
            for sign, shift in zip((-1, 1), cell_bounds[dim]):
                if P1[dim] == P2[dim] == center[dim]+shift:
                    return (dim, sign)
        raise ValueError('error!')

    cell_sizes = [b - a for a, b in cell_bounds]

    def comparator(P1, P2):
        lengths = 2*[None]

        def n_even_odd(n):
            even = sum(divmod(n, 2))
            return (even, n - even)

        for index, P in enumerate((P1, P2)):
            ds = get_dim_sign(P, center, cell_bounds)
            dim = ds[0]
            i = clockwise_ds.index(ds)
            n_even, n_odd = n_even_odd(i)
            # length0 = cell_sizes[0]*n_odd(i) + cell_sizes[1]*n_even(i)
            length0 = cell_sizes[0]*n_odd + cell_sizes[1]*n_even

            # length = 2*i + abs(P[1-dim]-cell_vertices[i][1-dim])
            lengths[index] = length0 + abs(P[1-dim]-cell_vertices[i][1-dim])

        return lengths[0] - lengths[1]

    side_pts = cell_vertices.copy()

    for index, subpts in enumerate(subpts_set):
        for side_pt in (subpts[0], subpts[-1]):
            side_pts.append(side_pt)

    side_pts = list(set(side_pts))
    loop = list(sorted(side_pts, key=cmp_to_key(comparator)))

    mid_pt = [0, 0]
    for P1, P2 in pairs(loop, last=True):
        dim, s = get_common_dim_sign(P1, P2)
        direction = [0, 0]
        # vogliamo gli indici spaziali invertiti
        direction[1-dim] = s

        mid_pt[dim] = P1[dim]
        mid_pt[1-dim] = (P1[1-dim]+P2[1-dim])/2
        if not is_inside_wall(mid_pt, pts):
            free_perimeter[tuple(direction)] += abs(P2[1-dim]-P1[1-dim])

    return free_perimeter


def is_inside_cell(P, center, cell_bounds):
    return all(P[dim] > center[dim] + cell_bounds[dim][0] and
               P[dim] < center[dim] + cell_bounds[dim][1]
               for dim in range(2))


def data_cells(pts, cell_bounds):
    checked_centers = []
    for P in pts:
        for shift in product(range(2), repeat=2):
            center = tuple(int(P[dim]+shift[dim]-P[dim] % 1)
                           for dim in range(2))
            if (center not in checked_centers and
               is_inside_cell(P, center, cell_bounds)):
                checked_centers.append(center)
                yield center


def get_subpts(center, pts, cell_bounds):

    def is_on_perimeter(P):
        return any(P[dim] == center[dim]+s and
                   P[1-dim] >= center[1-dim] + cell_bounds[1-dim][0] and
                   P[1-dim] <= center[1-dim] + cell_bounds[1-dim][1]
                   for dim in range(2) for s in cell_bounds[dim])

    def different_sides(P1, P2):
        if not (is_on_perimeter(P1) and is_on_perimeter(P2)):
            raise ValueError('Almeno uno dei due punti non è sul perimetro')
        return all(coord1 != coord2 for coord1, coord2 in zip(P1, P2))

    def following_points(index):
        n = 0
        while True:
            yield pts[(index+n) % len(pts)]
            n += 1

    subpts_set = []
    for index, (P1, P2) in enumerate(pairs(pts, last=True)):
        if is_on_perimeter(P1):
            if is_on_perimeter(P2) and different_sides(P1, P2):
                subpts_set.append([P1, P2])
            elif is_inside_cell(P2, center, cell_bounds):
                subpts = [P1, P2]

                # index corrisponde a P1, quindi dobbiamo mettere index+2
                for P in following_points(index+2):
                    subpts.append(P)
                    if is_on_perimeter(P):
                        break
                subpts_set.append(subpts)

    if len(subpts_set) > 2:
        # TODO: per ora lo lascierei così, ci sono problemi o
        # il reticolo diventa troppo fine lo togliamo..
        print(center)
        print(subpts_set)
        raise ValueError('Il reticolo non è abbastanza fine.')

    else:
        return subpts_set


def get_dim_sign(P, center, cell_bounds):
    ds_set = []
    for dim in range(2):
        for sign, shift in zip((-1, 1), cell_bounds[dim]):
            if P[dim] == center[dim]+shift:
                ds_set.append((dim, sign))

    if len(ds_set) > 1:
        for ds in clockwise_ds:
            if ds in ds_set:
                return ds
    elif len(ds_set) == 1:
        return ds_set[0]
    else:
        raise ValueError('error!')


def intersections_number(normal, P0, pts):
    # TODO: if possible, change or remove this part.
    global current_pts, c_pts
    if 'current_pts' not in globals() or pts != current_pts:
        c_pts = build_c_vecs(pts)
        current_pts = pts

    Vector2D = c_double*2
    c_normal = Vector2D(*normal)
    c_P0 = Vector2D(*P0)
    n_inters = (c_int*2)()
    n_pts = c_ulong(len(pts))
    inters_num_lib.intersections_number(n_inters, c_normal, c_P0, c_pts, n_pts)
    return list(n_inters)


def is_inside_wall(P0, pts):
    same_side, opp_side = intersections_number((0, 1), P0, pts)
    if same_side % 2 == opp_side % 2:
        return same_side % 2 != 0
    else:
        print('P0:', P0)
        raise ValueError('problem!')


def get_inward_normal_sign(normal, P0, pts):
    same_side, opp_side = intersections_number(normal, P0, pts)
    if same_side % 2 == opp_side % 2:
        return -1 if same_side % 2 == 0 else 1
    else:
        print('P0:', P0)
        print('normal:', normal)
        print('same_side, opp_side', same_side, opp_side)
        raise ValueError('problem!')
