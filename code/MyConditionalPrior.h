#ifndef _MyConditionalPrior_
#define _MyConditionalPrior_

#include "DNest4/code/RJObject/ConditionalPriors/ConditionalPrior.h"

class MyConditionalPrior:public DNest4::ConditionalPrior
{
	private:

		// Mean of amplitudes and widths
		double mu_loga, sigma_loga, mu_logq, sigma_logq;
		double pp;


		double perturb_hyperparameters(DNest4::RNG& rng);

	public:
		MyConditionalPrior();

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;

		static const int weight_parameter = 1;
};

#endif

