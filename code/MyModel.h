#ifndef _MyModel_
#define _MyModel_

#include "DNest4/code/RNG.h"
#include "MyConditionalPrior.h"
#include <ostream>
#include "Data.h"
#include "DNest4/code/RJObject/RJObject.h"
#include <vector>
 
class MyModel
{
	private:
		// Reference to the data
		static const Data& data;

		// A flat background level
		double background;

		// The Lorentzians
		DNest4::RJObject<MyConditionalPrior> dopplershift;
		// Extra white noise on the flux
		std::vector<double> noise_normals;
		double noise_sigma, noise_L;

		// Poisson mean
		std::vector<long double> mu;

		// Calculate mu from scratch:
		void calculate_mu();

	public: 
		// Constructor only gives size of params
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

		static double gaussian_cdf(double x, double x0, double gamma);
};

#endif

