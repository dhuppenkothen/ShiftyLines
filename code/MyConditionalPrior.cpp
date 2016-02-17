#include "MyConditionalPrior.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>

using namespace DNest4;
 
MyConditionalPrior::MyConditionalPrior()
{
}

void MyConditionalPrior::from_prior(RNG& rng)
{

	// Laplacian prior on the amplitude has parameters mu_loga and sigma_loga
	// mu_loga is uniformely distributed between -5 and 5
	mu_loga = (5. - (-5.))*rng.rand() + (-5.);

	// sigma_loga is uniformely distributed between 0 and 2
	sigma_loga = (2.-0.)*rng.rand() + 0.;


	// Laplacian prior on the log-q-factor has parameters mu_logq and sigma_logq
	// mu_logq is uniformely distributed between log(100) and log(1000)
	mu_logq = (log(1000.) - log(100.))*rng.rand() + log(100.);

	// sigma_logq is uniformely distributed between 0 and 0.3
	sigma_logq = (0.3 - 0.)*rng.rand() + 0.;

	// The parameter p deciding on the threshold for the signs has a Uniform distribution
	pp = (1. - 0.)*rng.rand() + 0.;

}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;

	int which = rng.rand_int(5);

	if(which == 0)
	{

		mu_loga += rng.randh()*(5.- (-5.));
                wrap(mu_loga, -5., 5.);
	}
	if(which == 1)
	{
                sigma_loga += rng.randh()*(2.- 0.);
                wrap(sigma_loga, 0., 2.);

	}
	if(which == 2)
	{
                mu_logq += rng.randh()*(log(1000.) - log(100.));
                wrap(mu_logq, log(100.), log(1000.));
	}
	if(which == 3)
	{
		sigma_logq += rng.randh()*(0.3 - 0.0);
		wrap(sigma_logq, 0.0, 0.3);
	}
	if(which == 4)
	{
		pp += rng.randh()*1.;
		wrap(pp, 0., 1.);
	}
	return logH;
}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{

	const int nlines = Data::get_instance().get_nlines();

	double loga, logq, sign;
	const double dshift = vec[0];

	if(dshift < -1. || dshift > 1.)
		return -1E300;

	double logprior = 0.;
	for (int i=0; i<nlines; i++)	
		{
			loga = vec[i+1];
			logq = vec[i+1+nlines];
			sign = vec[i+1+2*nlines];
	
			if(sign < 0.0 || sign > 1.)
				return -1E300;	
		
			logprior += -log(2.*sigma_loga) - std::abs(loga - mu_loga)/sigma_loga;
			logprior += -log(2.*sigma_logq) - std::abs(logq - mu_logq)/sigma_logq;

		}

	return logprior;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const

{

        const int nlines = Data::get_instance().get_nlines();

	vec[0] = -1. + (1. - (-1.))*vec[0]; // Doppler shift
	for (int i=0; i<nlines; i++)
		{
			if (vec[i+1] < 0.5)
				vec[i+1] = mu_loga + sigma_loga*log(2.*vec[i+1]);
			else
				vec[i+1] = mu_loga + sigma_loga*log(2. - 2.*vec[i+1]);

			if (vec[i+1+nlines] < 0.5)
				vec[i+1+nlines] = mu_logq + sigma_logq*log(2.*vec[i+1+nlines]);
			else
				vec[i+1+nlines] = mu_logq + sigma_logq*log(2. - 2.*vec[i+1+nlines]);

	
		}

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{

        const int nlines = Data::get_instance().get_nlines();

        vec[0] = 0.5*(vec[0] + 1.); // Doppler shift


        for (int i=0; i<nlines; i++)
                {
                        if (vec[i+1] < 0.5)
                        	vec[i+1] = 0.5*exp((vec[i+1] - mu_loga)/sigma_loga);
			else
				vec[i+1] = 1. - 0.5*exp((mu_loga - vec[i+1])/sigma_logq);
	
                        if (vec[i+1+nlines] < 0.5)
                                vec[i+1+nlines] = 0.5*exp((vec[i+1+nlines] - mu_loga)/sigma_loga);
                        else
                                vec[i+1+nlines] = 1. - 0.5*exp((mu_loga - vec[i+1+nlines])/sigma_logq);


                }

}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<mu_loga<<' '<<sigma_loga<<' '<<mu_logq<<' '<<sigma_logq<<' '<<pp<<' ';

}

