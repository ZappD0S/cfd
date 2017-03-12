#include <stdlib.h>
#include "spline_group.c"


#define IDX2D(i, j, nx) ((i)*(nx) + (j))

void swap(double** a, double** b) {
  double* tmp = *a;
  *a = *b;
  *b = tmp;
}


double unbounded_spline_group_eval(spline_group* spl_g, double y, double x) {
  if (y < 0 || y > (spl_g->ny - 1)) {
    return 0;
  } else if (x < 0) {
    return spline_group_eval(spl_g, y, 0);
  } else if (x > (spl_g->nx - 1)) {
    return spline_group_eval(spl_g, y, spl_g->nx - 1);
  } else {
    return spline_group_eval(spl_g, y, x);
  }
}


void advect_vels(spline_group** spl_g_set, double** uv, double dt)
{
  int dim, n;
  size_t i, j;
  double step_coeffs[] = {0.5, 0.5, 1};
  double k_coeffs[] = {6, 3, 3, 6};
  double k0[2], k[2], ds[2];
  double* p_k0 = k0;
  double* p_k = k;

  for (i = 0; i < spl_g_set[0]->ny; i++) {
    for (j = 1; j < spl_g_set[0]->nx - 1; j++) {

      for (dim = 0; dim < 2; dim++) {
        p_k[dim] = unbounded_spline_group_eval(spl_g_set[dim], i, j);
        ds[dim] = p_k[dim]/k_coeffs[0];
      }

      for (n = 1; n < 4; n++) {
        swap(&p_k, &p_k0);

        for (dim = 0; dim < 2; dim++) {
          p_k[dim] = unbounded_spline_group_eval(
            spl_g_set[dim],
            i - dt*step_coeffs[n-1]*p_k0[1],
            j - dt*step_coeffs[n-1]*p_k0[0]);
          ds[dim] += p_k[dim]/k_coeffs[n];
          }
        }

      for (dim = 0; dim < 2; dim++) {
        uv[dim][IDX2D(i, j, spl_g_set[0]->nx)] = unbounded_spline_group_eval(
          spl_g_set[dim], i - dt*ds[1], j - dt*ds[0]);
      }
    }
  }
}
