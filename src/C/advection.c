#include <stdlib.h>
#include "splines.c"
// #include "spline_group.c"

// gcc -fPIC -O2 advection.c -shared -o ../../bin/advection.so

#define swap(a, b) do {   \
  double* tmp = *a;       \
  *a = *b;                \
  *b = tmp;               \
} while(0)


#define IDX2D(spl, i, j) (spl->nx * (i) + (j))


void advect_vels(Spline2DGroup* spl_g, double** uv, double dt)
{
  int dim, n;
  size_t i, j;
  const double step_coeffs[] = {0.5, 0.5, 1};
  const double k_coeffs[] = {6, 3, 3, 6};
  double k0[2], k[2], ds[2];
  double* k_ptr = k;
  double* k0_ptr = k0;

  for (i = 0; i < spl_g->ny; ++i) {
    for (j = 1; j < spl_g->nx - 1; ++j) {
      Spline2DGroup_unbounded_eval(spl_g, j, i, k);

      for (dim = 0; dim < 2; ++dim)
        ds[dim] = k[dim]/k_coeffs[0];

      for (n = 1; n < 4; ++n) {
        swap(&k_ptr, &k0_ptr);
        Spline2DGroup_unbounded_eval(
          spl_g,
          j - dt*step_coeffs[n-1]*k0[0],
          i - dt*step_coeffs[n-1]*k0[1],
          k);

        for (dim = 0; dim < 2; ++dim)
          ds[dim] += k[dim]/k_coeffs[n];
      }

      for (dim = 0; dim < 2; ++dim)
        Spline2D_unbounded_eval(
          spl_g->splines[dim],
          j - dt*ds[0],
          i - dt*ds[1],
          uv[dim] + IDX2D(spl_g, i, j));
    }
  }
}


// double unbounded_spline_group_eval(spline_group* spl_g, double y, double x) {
//   if (y < 0 || y > (spl_g->ny - 1)) {
//     return 0;
//   } else if (x < 0) {
//     return spline_group_eval(spl_g, y, 0);
//   } else if (x > (spl_g->nx - 1)) {
//     return spline_group_eval(spl_g, y, spl_g->nx - 1);
//   } else {
//     return spline_group_eval(spl_g, y, x);
//   }
// }
//
// void advect_vels(spline_group** spl_g_set, double** uv, double dt)
// {
//   int dim, n;
//   size_t i, j;
//   double step_coeffs[] = {0.5, 0.5, 1};
//   double k_coeffs[] = {6, 3, 3, 6};
//   double k[2], k0[2], ds[2];
//   double* k_ptr = k;
//   double* k0_ptr = k0;
//
//   for (i = 0; i < spl_g_set[0]->ny; i++) {
//     for (j = 1; j < spl_g_set[0]->nx - 1; j++) {
//
//       for (dim = 0; dim < 2; dim++) {
//         k[dim] = unbounded_spline_group_eval(spl_g_set[dim], i, j);
//         ds[dim] = k[dim]/k_coeffs[0];
//       }
//
//       for (n = 1; n < 4; n++) {
//         swap(&k_ptr, &k0_ptr);
//
//         for (dim = 0; dim < 2; dim++) {
//           k[dim] = unbounded_spline_group_eval(
//             spl_g_set[dim],
//             i - dt*step_coeffs[n-1]*k0[1],
//             j - dt*step_coeffs[n-1]*k0[0]);
//           ds[dim] += k[dim]/k_coeffs[n];
//           }
//         }
//
//       for (dim = 0; dim < 2; dim++) {
//         uv[dim][IDX2D(spl_g_set[0], i, j)] = unbounded_spline_group_eval(
//           spl_g_set[dim], i - dt*ds[1], j - dt*ds[0]);
//       }
//     }
//   }
// }
