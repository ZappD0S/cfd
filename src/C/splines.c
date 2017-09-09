#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// gcc -fPIC -O2 splines.c -shared -o ../../bin/splines.so

#define IDX2D(spl, i, j) (spl->nx * (i) + (j))

#define cubic_eval(p, x)																		\
((p)[1]																											\
	+ (x)*((p)[2] - (p)[0]																		\
		+ (x)*(2*(p)[0] - 5*(p)[1] + 4*(p)[2] - (p)[3]					\
			+ (x)*(3*((p)[1] - (p)[2]) + (p)[3] - (p)[0])))/2)


typedef struct Spline2D {
	size_t ny;
	size_t nx;
	double* data;
	void (*eval)(struct Spline2D* spl, double* point,  double* z);
} Spline2D;


void bilinear_eval(Spline2D* spl, double* point,  double* z)
{
	double x = point[0], y = point[1];
	int xi = (int) x, yi = (int) y;

	if (x == xi && y == yi) {
		*z = spl->data[IDX2D(spl, yi, xi)];
		return;
	}
	int i, j;
	double xfrac, yfrac;
	double dx = x - xi, dy = y - yi;
	*z = 0;
	for (i = 0; i < 2; ++i) {
  	for (j = 0; j < 2; ++j) {
    	xfrac = (j) ? dx : 1. - dx;
    	yfrac = (i) ? dy : 1. - dy;
    	*z += xfrac * yfrac * spl->data[IDX2D(spl, yi + i, xi + j)];
  	}
	}
}


void bicubic_eval(Spline2D* spl, double* point,  double* z)
{
	double x = point[0], y = point[1];
	int xi = (int) x, yi = (int) y;

	if (x == xi && y == yi) {
		*z = spl->data[IDX2D(spl, yi, xi)];
		return;
	} else if (xi < 2 || xi > (spl->nx - 1 - 3) ||
						 yi < 2 || yi > (spl->ny - 1 - 3)) {
		bilinear_eval(spl, point, z);
		return;
	}

	xi -= 1;
	yi -= 1;

	assert((x - xi - 1) >= 0 && (x - xi - 1) <= 1);
	assert((y - yi - 1) >= 0 && (y - yi - 1) <= 1);

	double arr[4];
	int i;
	for (i = 0; i < 4; ++i)
		arr[i] = cubic_eval(spl->data + IDX2D(spl, yi + i, xi), x - xi - 1);
	*z = cubic_eval(arr, y - yi - 1);
}


Spline2D* Spline2D_alloc(size_t* shape, int type)
{
	Spline2D* spl = (Spline2D*) malloc(sizeof(Spline2D));
	spl->ny = shape[0];
	spl->nx = shape[1];
  spl->data = (double*) malloc(spl->ny * spl->nx * sizeof(double));
	switch (type) {
		case 1:
			spl->eval = &bilinear_eval;
			break;
		case 2:
			spl->eval = &bicubic_eval;
			break;
	}
	return spl;
}


void Spline2D_set(Spline2D* spl, double* data)
{
	int i, j;
	for (i = 0; i < spl->ny; ++i) {
		for (j = 0; j < spl->nx; ++j) {
			spl->data[spl->nx * i + j] = data[spl->nx * i + j];
		}
	}
}


void Spline2D_eval(Spline2D* spl, double* point,  double* z)
{
	spl->eval(spl, point, z);
}


void Spline2D_unbounded_eval(Spline2D* spl, double x, double y, double* z)
{
	if (y < 0 || y > (spl->ny - 1)) {
		*z = 0;
	} else {
		double point[2] = {x, y};
		if (x < 0) {
			point[0] = 0;
		} else if (x > (spl->nx - 1)) {
			point[0] = spl->nx - 1;
		}
		Spline2D_eval(spl, point, z);
	}
}


void Spline2D_array_eval(
	Spline2D* spl, double** points,  double* z_arr, size_t size)
{
	size_t i;
	for (i = 0; i < size; ++i)
		Spline2D_eval(spl, points[i] , z_arr + i);
}


void Spline2D_free(Spline2D* spl)
{
	free(spl->data);
	free(spl);
}


typedef struct {
	Spline2D** splines;
	int splines_number;
	size_t nx;
	size_t ny;
} Spline2DGroup;


Spline2DGroup* Spline2DGroup_alloc(size_t* shape, int type, int n)
{
	Spline2DGroup* spl_g = (Spline2DGroup*) malloc(sizeof(Spline2DGroup));
	spl_g->splines = (Spline2D**) malloc(n * sizeof(Spline2D*));
	int i;
	for (i = 0; i < n; ++i)
		spl_g->splines[i] = Spline2D_alloc(shape, type);

	spl_g->splines_number = n;
	spl_g->ny = shape[0];
	spl_g->nx = shape[1];
	return spl_g;
}


void Spline2DGroup_set(Spline2DGroup* spl_g, double** data_arr)
{
  size_t i;
  for (i = 0; i < spl_g->splines_number; ++i)
  	Spline2D_set(spl_g->splines[i], data_arr[i]);
}


void Spline2DGroup_eval(
	Spline2DGroup* spl_g, double* point, double* z_arr)
{
  size_t i;
  for (i = 0; i < spl_g->splines_number; ++i)
		Spline2D_eval(spl_g->splines[i], point, z_arr + i);
}


void Spline2DGroup_unbounded_eval(
	Spline2DGroup* spl_g, double x, double y, double* z_arr)
{
	size_t i;
	for (i = 0; i < spl_g->splines_number; ++i)
		Spline2D_unbounded_eval(spl_g->splines[i], x, y, z_arr + i);
}


void Spline2DGroup_array_eval(
	Spline2DGroup* spl_g, double** points,
	double** z_arrs, size_t size)
{
  size_t i;
  for (i = 0; i < spl_g->splines_number; ++i)
  	Spline2D_array_eval(spl_g->splines[i], points, z_arrs[i], size);
}


void Spline2DGroup_free(Spline2DGroup* spl_g) {
    size_t i;
    for (i = 0; i < spl_g->splines_number; ++i)
        Spline2D_free(spl_g->splines[i]);

    free(spl_g->splines);
    free(spl_g);
}
