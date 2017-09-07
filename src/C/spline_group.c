#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>


//#define IDX2D(i, j, nx) ((i)*(nx) + (j))
#define IDX2D(spl, i, j) (spl->nx * (i) + (j))

// gcc -shared -o spline_group.so -fPIC spline_group.c `gsl-config --cflags --libs`


typedef struct spline_group {
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  double *za;
  double *xa;
  double *ya;
  size_t nx;
  size_t ny;
} spline_group;


spline_group* alloc_spline_group(size_t ny, size_t nx) {
  spline_group* spl_g = (spline_group*) calloc(1, sizeof(spline_group));
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  spl_g->spline = gsl_spline2d_alloc(T, ny, nx);
  size_t n;
  spl_g->xa = malloc(nx * sizeof(double));
  for (n = 0; n < nx; n++) {
    spl_g->xa[n] = (double) n;
  }

  spl_g->ya = malloc(ny * sizeof(double));
  for (n = 0; n < ny; n++) {
    spl_g->ya[n] = (double) n;
  }

  spl_g->za = malloc(nx * ny * sizeof(double));

  spl_g->xacc = gsl_interp_accel_alloc();
  spl_g->yacc = gsl_interp_accel_alloc();

  spl_g->nx = nx;
  spl_g->ny = ny;

  return spl_g;
}


void init_spline_group(spline_group* spl_g, double* u) {
  size_t i, j;
  for (i = 0; i < spl_g->ny; i++) {
    for (j = 0; j < spl_g->nx; j++) {
      gsl_spline2d_set(spl_g->spline, spl_g->za, i, j,
                       u[IDX2D(spl_g, i, j)]);
    }
  }
  gsl_spline2d_init(spl_g->spline, spl_g->ya, spl_g->xa, spl_g->za,
                    spl_g->ny, spl_g->nx);
}


double spline_group_eval(spline_group* spl_g, double y, double x) {
  return gsl_spline2d_eval(spl_g->spline, y, x, spl_g->yacc, spl_g->xacc);
}


void spline_group_array_eval(spline_group* spl_g, double* y, double** x,
                             size_t size)
{
  size_t n;
  for (n = 0; n < size; n++) {
    y[n] = spline_group_eval(spl_g, x[n][1], x[n][0]);
    // printf("res: %f\n", y[n]);
  }
}


void free_spline_group(spline_group* spl_g) {
  gsl_spline2d_free(spl_g->spline);
  gsl_interp_accel_free(spl_g->xacc);
  gsl_interp_accel_free(spl_g->yacc);
  free(spl_g->xa);
  free(spl_g->ya);
  free(spl_g->za);
  free(spl_g);
}
