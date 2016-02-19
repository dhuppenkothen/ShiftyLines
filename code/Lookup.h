#ifndef ShiftyLines_Lookup
#define ShiftyLines_Lookup

#include <vector>

/*
* Lookup tables for speeding things up
* Singleton pattern
*/

class Lookup
{
	private:
		const std::size_t num;
		const double x_min;
		const double x_max;
		const double dx;
		const double one_over_dx;

		std::vector<double> evaluations;

		Lookup(std::size_t num, double x_min, double x_max);
		Lookup(const Lookup& other);

		static Lookup instance;

	public:
		double evaluate_erf(double x);
		static double erf(double x);

};

#endif

