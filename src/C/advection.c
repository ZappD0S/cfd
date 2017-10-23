#include <stdlib.h>
#include "splines.c"

// gcc -fPIC -O2 advection.c -shared -o ../../bin/advection.so

void swap(double** a, double** b) {
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}


void advect_vels(Spline2DGroup* spl_g, double** uv, double dt)
{
  int dim, n;
  size_t i, j;
  const double step_coeffs[] = {0.5, 0.5, 1};
  const double k_coeffs[] = {6, 3, 3, 6};
  double k0_arr[2], k_arr[2], ds[2];
  double* k = k_arr;
  double* k0 = k0_arr;

  for (i = 1; i < spl_g->ny - 1; ++i) {
    for (j = 1; j < spl_g->nx - 1; ++j) {
      Spline2DGroup_unbounded_eval(spl_g, j, i, k);

      for (dim = 0; dim < 2; ++dim)
        ds[dim] = k[dim]/k_coeffs[0];

      for (n = 1; n < 4; ++n) {
        swap(&k, &k0);
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
