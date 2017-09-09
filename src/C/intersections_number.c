#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// gcc -fPIC -O2 intersections_number.c -shared -o ../../bin/intersections_number.so

#define dot2D(v, w) ((v[0])*(w[0]) + (v[1])*(w[1]))
#define cross2D(v, w) ((v[0])*(w[1]) - (v[1])*(w[0]))
#define opp_signs(a, b) ((a > 0 && b < 0) || (a < 0 && b > 0))


void intersections_number(int *n_inters, double *normal, double *origin, double **pts, size_t n_pts)
{
  double P0[2], P1[2], P2[2];
  size_t i;
  int j, a, b, on_border = 0;

  for (i = 0; i < n_pts; ++i) {
    for (j = 0; j < 2; ++j) {
      P0[j] = pts[(i - 1) % n_pts][j] - origin[j];
      P1[j] = pts[i][j] - origin[j];
      P2[j] = pts[(i + 1) % n_pts][j] - origin[j];
    }

    if (cross2D(P0, P1) == 0 || cross2D(P1, P2) == 0) {
      on_border = 1;
    } else {
      a = opp_signs(cross2D(P1, normal), cross2D(P2, normal));
      b = (cross2D(P1, normal) == 0 &&
           opp_signs(cross2D(P0, normal), cross2D(P2, normal)));
      if (a || b) {
        if ((a && (dot2D(P1, normal) > 0 || dot2D(P2, normal) > 0)) ||
            (b && dot2D(P1, normal) > 0)) {
          n_inters[0] += 1;
        } else {
          n_inters[1] += 1;
        }
      }
    }
  }
  if (on_border && (n_inters[0] || n_inters[1])) {
    for (j = 0; j < 2; j++) {
      if (n_inters[j] % 2 == 0) {
        n_inters[j] += 1;
        break;
      }
    }
  }
}


// void intersections_number(int *n_inters, double *normal, double *origin, double **pts, size_t n_pts)
// {
//   double P0[2], P1[2], P2[2];
//   size_t i, j;
//   int on_border = 0;
//
//   for (i = 0; i < n_pts; i++) {
//     for (j = 0; j < 2; j++) {
//       P0[j] = pts[(i - 1) % n_pts][j] - origin[j];
//       P1[j] = pts[i][j] - origin[j];
//       P2[j] = pts[(i + 1) % n_pts][j] - origin[j];
//     }
//     if (cross2D(P0, P1) == 0 || cross2D(P1, P2) == 0 || P1) == 0 {
//       on_border = 1;
//     } else if (opp_signs(cross2D(P1, normal), cross2D(P2, normal))) {
//       if (dot2D(P1, normal) > 0 || dot2D(P2, normal) > 0) {
//         n_inters[0] += 1;
//       } else {
//         n_inters[1] += 1;
//       }
//     } else if (cross2D(P1, normal) == 0 &&
//                opp_signs(cross2D(P0, normal), cross2D(P2, normal))) {
//       if (dot2D(P0, normal) > 0 ||
//           dot2D(P1, normal) > 0 ||
//           dot2D(P2, normal) > 0) {
//         n_inters[0] += 1;
//       } else {
//         n_inters[1] += 1;
//       }
//     }
//   }
//   if (on_border && (n_inters[0] || n_inters[1])) {
//     for (j = 0; j < 2; j++) {
//       if (n_inters[j] % 2 == 0) {
//         n_inters[j] += 1;
//         break;
//       }
//     }
//   }
// }
