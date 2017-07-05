#include "MyConditionalPrior.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>

using namespace DNest4;
 

// These are the lower and upper bounds on the uniform 
// prior on the Doppler shift. Change them here if you'd 
// like the prior to change!
const double MyConditionalPrior::dmin = -0.01;
const double MyConditionalPrior::dmax = 0.01;

 
MyConditionalPrior::MyConditionalPrior()
{
}

void MyConditionalPrior::from_prior(RNG& rng)
{

	// NOTE: We want to change the prior on the mean of the Laplacian 
	// prior on the amplitude to a Cauchy distribution, because our 
	// amplitudes tend to be smaller than -5


	// Laplacian prior on the amplitude has parameters mu_loga and sigma_loga
	// mu_loga is uniformely distributed between -5 and 5
//	mu_loga = (5. - (-5.))*rng.rand() + (-5.);
        mu_loga = tan(M_PI*(0.97*rng.rand() - 0.485));

	// sigma_loga is uniformely distributed between 0 and 2
	sigma_loga =  (4.-0.)*rng.rand() + 0.;


	// Laplacian prior on the log-q-factor has parameters mu_logq and sigma_logq
	// mu_logq is uniformely distributed between log(100) and log(1000)
	//mu_logq = (log(1000.) - log(300.))*rng.rand() + log(300.);
	double min_df = Data::get_instance().get_min_df(); 
	double f_range = Data::get_instance().get_f_range() ;
	mu_logwidth = (log(f_range/10.) - log(min_df/5.0))*rng.rand() + log(min_df/5.0);
  
	// sigma_logq is uniformely distributed between 0 and 0.3
	sigma_logwidth =  (1.0 - 0.)*rng.rand() + 0.;

	// The parameter p deciding on the threshold for the signs has a Uniform distribution
	pp = (1.0 - 0.)*rng.rand() + 0.;

        // The parameter p deciding on the threshold for the signs has a Uniform distribution
        pp_presence = (1.0 - 0.)*rng.rand() + 0.;
        
}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;
        double min_df = Data::get_instance().get_min_df();
        double f_range = Data::get_instance().get_f_range() ;

	int which = rng.rand_int(6);

	if(which == 0)
	{
//		mu_loga += rng.randh()*(5.- (-5.));
//                wrap(mu_loga, -5., 5.);
                mu_loga = (atan(mu_loga)/M_PI + 0.485)/0.97;
                mu_loga += pow(10., 1.5 - 6.*rng.rand())*rng.randn();
                mu_loga = mod(mu_loga, 1.);
                mu_loga = tan(M_PI*(0.97*mu_loga - 0.485));
      


	}
	if(which == 1)
	{
                sigma_loga += rng.randh()*(4.- 0.);
                wrap(sigma_loga, 0., 4.);
	}
	if(which == 2)
	{
                mu_logwidth += rng.randh()*(log(f_range/10.) - log(min_df/5.0));
                wrap(mu_logwidth, log(min_df/5.0), log(f_range/10.));
	}
	if(which == 3)
	{
		sigma_logwidth += rng.randh()*(1.0 - 0.0);
		wrap(sigma_logwidth, 0.0, 1.0);
	}
	if(which == 4)
	{
		pp += rng.randh()*1.0;
		wrap(pp, 0., 1.0);
	}
        if(which == 4)
        {
                pp_presence += rng.randh()*1.0;
                wrap(pp_presence, 0., 1.0);
        }

	return logH;
}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{

	const int nlines = Data::get_instance().get_nlines();

	double loga, logwidth, the_sign, presence;
	const double dshift = vec[0];

	if(dshift < dmin || dshift > dmax)
		return -1E300;

	double logprior = 0.;

//	logprior += -log(0.01) - log(1. + pow(dshift, 2)/pow(0.01, 2));

	for (int i=0; i<nlines; i++)	
		{
			loga = vec[i+1];
			logwidth = vec[i+1+nlines];
			the_sign = vec[i+1+2*nlines];
			presence = vec[i+1+3*nlines];
	
			if(the_sign < 0.0 || the_sign > 1.)
				return -1E300;	

                        if(presence < 0.0 || presence > 1.)
                                return -1E300; 

			//if(exp(loga) > background)
			//	return -1E300;
		
			logprior += -log(2.*sigma_loga) - std::abs(loga - mu_loga)/sigma_loga;
			logprior += -log(2.*sigma_logwidth) - std::abs(logwidth - mu_logwidth)/sigma_logwidth;

		}

	return logprior;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const

{

        const int nlines = Data::get_instance().get_nlines();

	vec[0] = dmin + (dmax - dmin)*vec[0]; // Doppler shift
//	vec[0] = 0.01*tan(M_PI*(vec[0] - 0.5));

	
	for (int i=0; i<nlines; i++)
		{
//			if (vec[i+1] < 0.5)
//				vec[i+1] = mu_loga + sigma_loga*log(2.*vec[i+1]);
//			else
//				vec[i+1] = mu_loga - sigma_loga*log(2. - 2.*vec[i+1]);

				vec[i+1] = laplacian_cdf_inverse(vec[i+1], mu_loga, sigma_loga);
	
//			if (vec[i+1+nlines] < 0.5)
//				vec[i+1+nlines] = mu_logq + sigma_logq*log(2.*vec[i+1+nlines]);
//			else
//				vec[i+1+nlines] = mu_logq - sigma_logq*log(2. - 2.*vec[i+1+nlines]);

				vec[i+1+nlines] = laplacian_cdf_inverse(vec[i+1+nlines], mu_logwidth, sigma_logwidth);

	
		}

}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{

        const int nlines = Data::get_instance().get_nlines();
 
	vec[0] = (vec[0] - dmin)/(dmax - dmin);
//	vec[0] = (1./M_PI)*atan(vec[0]/0.01) + 0.5;


        for (int i=0; i<nlines; i++)
                {
			vec[i+1] = laplacian_cdf(vec[i+1], mu_loga, sigma_loga);
			vec[i+1+nlines] = laplacian_cdf(vec[i+1+nlines], mu_logwidth, sigma_logwidth);
//                        if (vec[i+1] < mu_loga)
//                        	vec[i+1] = 0.5*exp((vec[i+1] - mu_loga)/sigma_loga);
//			else
//				vec[i+1] = 1. - 0.5*exp((mu_loga - vec[i+1])/sigma_loga);
	
//                        if (vec[i+1+nlines] < mu_logq)
//                                vec[i+1+nlines] = 0.5*exp((vec[i+1+nlines] - mu_logq)/sigma_logq);
//                        else
//                                vec[i+1+nlines] = 1. - 0.5*exp((mu_logq - vec[i+1+nlines])/sigma_logq);


                }

}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<mu_loga<<' '<<sigma_loga<<' '<<mu_logwidth<<' '<<sigma_logwidth<<' '<<pp<<' '<<pp_presence<<' ';
}


// Laplacian cdf stuff
int MyConditionalPrior::sign(double x)
{
	if(x == 0.)
		return 0;
	if(x > 0.)
		return 1;
	return -1;
}

double MyConditionalPrior::laplacian_cdf(double x, double center, double width)
{
	assert(width > 0.);
	return 0.5 + 0.5*sign(x - center)*(1. - exp(-std::abs(x - center)/width));
}

double MyConditionalPrior::laplacian_cdf_inverse(double x, double center,
														double width)
{
	assert(width > 0.);
	assert(x >= 0. && x <= 1.);

	return center - width*sign(x - 0.5)*log(1. - 2.*std::abs(x - 0.5));
}

