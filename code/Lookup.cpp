#include "Lookup.h"
#include <cmath>
#include <cassert>

using namespace std;

Lookup Lookup::instance(1001, -10., 10.);

Lookup::Lookup(size_t num, double x_min, double x_max)
:num(num)
,x_min(x_min)
,x_max(x_max)
,dx((x_max - x_min)/(num - 1))
,one_over_dx(1./dx)
,evaluations(num)
{
	assert(x_max > x_min);

	for(size_t i=0; i<num; ++i)
		evaluations[i] = std::erf(x_min + i*dx);
}

double Lookup::evaluate_erf(double x)
{
	if(x < Lookup::x_min)
		return -1.;
	if(x > Lookup::x_max)
		return 1.;

	size_t index = static_cast<size_t>((x - x_min)*one_over_dx);
	double lambda = (x - (x_min + index*dx))*one_over_dx;

	return (1. - lambda)*evaluations[index] + lambda*evaluations[index+1];
}

double Lookup::erf(double x)
{
	return instance.evaluate_erf(x);
}

