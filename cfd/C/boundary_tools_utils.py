from ctypes import c_double, c_ulong, c_int, POINTER, pointer, CDLL
import os

pwd = os.path.dirname(__file__)
lib = CDLL(os.path.join(pwd, 'intersections_number.so'))


def build_c_pts(pts):
    Vector2D = c_double*2
    c_arr = (POINTER(Vector2D)*len(pts))()
    for n, P in enumerate(pts):
        c_arr[n] = pointer(Vector2D(*P))
    return pointer(c_arr)


def intersections_number(normal, P0, pts):
    global current_pts, c_pts
    if 'current_pts' not in globals() or pts != current_pts:
        c_pts = build_c_pts(pts)
        current_pts = pts

    Vector2D = c_double*2
    c_normal = Vector2D(*normal)
    c_P0 = Vector2D(*P0)
    n_inters = (c_int*2)()
    n_pts = c_ulong(len(pts))
    lib.intersections_number(n_inters, c_normal, c_P0, c_pts, n_pts)
    return list(n_inters)
