#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

// gcc -fPIC -O2 intersections_number.c -shared -o ../../bin/intersections_number.so

#define dot2D(v, w) ((v)[0]*(w)[0] + (v)[1]*(w)[1])
#define cross2D(v, w) ((v)[0]*(w)[1] - (v)[1]*(w)[0])
#define opp_signs(a, b) (((a) > 0 && (b) < 0) || ((a) < 0 && (b) > 0))
#define parallels(v, w) (cross2D(v, w) == 0)


void
intersections_number(
  int *n_inters, double *normal, double *origin, double **pts, size_t n_pts)
{
  double P0[2], P1[2], P2[2];
  double cross0, cross1, cross2;
  size_t i;
  int j, on_border = 0;

  for (i = 0; i < n_pts; ++i) {
    for (j = 0; j < 2; ++j) {
      P0[j] = pts[(i - 1) % n_pts][j] - origin[j];
      P1[j] = pts[i][j] - origin[j];
      P2[j] = pts[(i + 1) % n_pts][j] - origin[j];
    }

    if (cross2D(P1, P2) == 0) {
      on_border = 1;
      continue;
    }

    cross0 = cross2D(P0, normal);
    cross1 = cross2D(P1, normal);
    cross2 = cross2D(P2, normal);
    if (cross1 == 0) {
      assert(cross2D(P0, P2) != 0);
    }

    if (opp_signs(cross1, cross2)) {
      if (dot2D(P1, normal) > 0 || dot2D(P2, normal) > 0) {
        n_inters[0] += 1;
      } else {
        n_inters[1] += 1;
      }
    } else if (cross2D(P0, P1) != 0 && cross1 == 0 && opp_signs(cross0, cross2)) {
      if (dot2D(P0, normal) > 0 || dot2D(P2, normal) > 0) {
        n_inters[0] += 1;
      } else {
        n_inters[1] += 1;
      }
    }
  }

  if (on_border && (n_inters[0] || n_inters[1])) {
    for (j = 0; j < 2; ++j) {
      if ((n_inters[j] % 2) == 0) {
        n_inters[j] += 1;
        break;
      }
    }
  }
}


void
is_inside(int *res, double *origin, double **pts, size_t n_pts)
{
  int n_inters[2] = {0, 0};
  double normal[2] = {0, 1};
  intersections_number(n_inters, normal, origin, pts, n_pts);
  int same_side = n_inters[0], opp_side = n_inters[1];
  if ((same_side % 2) == (opp_side % 2)) {
    *res = (same_side % 2) == 1;
  } else {
    printf("same_side: %d, opp_side: %d\n", same_side, opp_side);
    printf("[%f, %f]\n", origin[0], origin[1]);
    assert(0);
  }
}


void is_inside2(int* res, double* P, double** pts, size_t n_pts)
{
  double *P1, *P2;
  int i, j;
  for (i = 0, j = n_pts-1; i < n_pts; j = i++) {
    P1 = pts[i];
    P2 = pts[j];

    if (((P1[1] > P[1]) != (P2[1] > P[1])) &&
        (P[0] < (P2[0] - P1[0]) * (P[1] - P1[1]) / (P2[1] - P1[1]) + P1[0]))
        *res = !(*res);
  }
}


// int
// ray_intersection(double* normal, double* origin, double** pts, size_t n_pts)
// {
//   double P0[2], P1[2], P2[2];
//   double theta0, theta1, theta2, theta3;
//   double cross0, cross1, cross2;
//   size_t i;
//   int j, net_inters = 0;
//
//   printf("starting\n");
//   for (i = 0; i < n_pts; ++i) {
//     for (j = 0; j < 2; ++j) {
//       P0[j] = pts[(i - 1) % n_pts][j] - origin[j];
//       P1[j] = pts[i][j] - origin[j];
//       P2[j] = pts[(i + 1) % n_pts][j] - origin[j];
//     }
//
//     if (parallels(P1, normal)) {
//       for (j = 0; j < 2; ++j)
//         P1[j] = pts[(i - 1) % n_pts][j] - origin[j];
//     } else if (parallels(P2, normal)) {
//       continue;
//     }
//
//     theta1 = atan2(P1[1], P1[0]);
//     theta2 = atan2(normal[1], normal[0]);
//     theta3 = atan2(P2[1], P2[0]);
//     theta0 = fmin(fmin(theta1, theta2), theta3);
//     theta1 -= theta0;
//     theta2 -= theta0;
//     theta3 -= theta0;
//     assert(theta1 != theta2 && theta2 != theta3);
//
//     if (theta1 < theta2 && theta2 < theta3) {
//       net_inters += 1;
//       printf("net_inters: %d\n", net_inters);
//     } else if (theta3 < theta2 && theta2 < theta1) {
//     // } else if (theta1 > theta2 && theta2 > theta3) {
//       net_inters -= 1;
//       printf("net_inters: %d\n", net_inters);
//     }
//   }
//   printf("end\n\n");
//   return net_inters;
// }
//
// void
// line_intersections(
//   int* net_inters, double* normal, double* origin, double** pts, size_t n_pts)
// {
//   int dir, sign, i;
//   double dir_normal[2];
//
//   for (dir = 0; dir < 2; ++dir) {
//     sign = -2*dir + 1;
//     for (i = 0; i < 2; ++i)
//       dir_normal[i] = sign*normal[i];
//
//     net_inters[dir] = ray_intersection(dir_normal, origin, pts, n_pts);
//   }
// }
