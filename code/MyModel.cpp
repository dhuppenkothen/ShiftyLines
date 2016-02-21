#include "MyModel.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include "Lookup.h"
#include <cmath>

using namespace std;
using namespace DNest4;

const Data& MyModel::data = Data::get_instance();
const int nlines = Data::get_instance().get_nlines();
//const std::vector<double> line_pos = Data::get_instance().get_line_pos();

MyModel::MyModel()
//:dopplershift(3*nlines+1, 5, false, MyConditionalPrior())
:dopplershift(3*nlines+1, 3, false, MyConditionalPrior())
,noise_normals(data.get_f_left().size())
,mu(data.get_f_left().size())
{
}

double MyModel::gaussian_cdf(double x, double x0, double gamma)
{
	return 0.5*(1. + Lookup::erf((x-x0)/(gamma*sqrt(2.))));
}

void MyModel::calculate_mu()
{ 
	// define line positions
	const vector<double>& line_pos = data.get_line_pos();
	const int nlines = data.get_nlines();

	// declare shifted line positions
	vector<double> line_pos_shifted;

	// get left and right boundaries of the wavelength bins
        const vector<double>& f_left = data.get_f_left();
        const vector<double>& f_right = data.get_f_right();
//        const vector<double>& f_mid = data.get_f_mid();

	// assign constant background to model
	mu.assign(mu.size(), background); //old version
 
	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();

	vector<double> amplitude, logq, sign, width;
	double dshift;

	vector<double> mu_temp(mu.size(), 0.0);

        for(size_t j=0; j<dopplershiftcomponents.size(); j++)
        {
	        // assign 0 to all shifted line positions
	        line_pos_shifted.assign(line_pos.size(), 0.0);

		amplitude.assign(nlines, 0.);
		logq.assign(nlines, 0.);
		sign.assign(nlines, 0.);
		width.assign(nlines, 0.);

		dshift = dopplershiftcomponents[j][0];
	
		for (int i=0; i<nlines; i++)
			{
				line_pos_shifted[i] = line_pos[i]*(1. + dshift);
				amplitude[i] = exp(dopplershiftcomponents[j][i+1]);
				logq[i] = dopplershiftcomponents[j][i+1+nlines];
				sign[i] = dopplershiftcomponents[j][i+1+2*nlines];		
				width[i] = line_pos_shifted[i]/exp(logq[i]);
			}
	
		
		int s=0;	
                for(size_t i=0; i<mu.size(); i++)
                {
			for (int k=0; k<nlines; k++)
				{
			// Integral over the Lorentzian distribution
					if (sign[k] < pp) 
						s = -1;
					else 
						s = 1;
 
				//	mu_temp[i] += s*amplitude[k]/(width[k]*sqrt(2.*M_PI))*exp(-pow(f_mid[i]-line_pos_shifted[k],2)/(2.*pow(width[k],2)));
					mu_temp[i] += s*amplitude[k]*(gaussian_cdf(f_right[i], line_pos_shifted[k], width[k])
								- gaussian_cdf(f_left[i], line_pos_shifted[k], width[k]));

				}

                } 

	}

        // Compute the OU process
        vector<double> y(mu.size());
        double alpha = exp(-1./noise_L);
        // y[i+1] = a*y[i] + sigma*n[i]
        // S^2 = a^2 S^2 + sigma^2
        // S^2 = sigma^2/(1 - a^2)
        // S = sigma/sqrt(1 - a^2)

        for(size_t i=0; i<mu.size(); i++)
        {
                if(i==0)
                        y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
                else
                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
                mu[i] *= exp(y[i]);
        }


        
	for (size_t i=0; i<mu.size(); i++)
		{
			mu[i]  += mu_temp[i];//exp(mu_temp[i]);
			if (mu[i] < 0.0)
				mu[i] = 0.0;
		}
}

void MyModel::from_prior(RNG& rng)
{
	background = tan(M_PI*(0.97*rng.rand() - 0.485));
	background = exp(background);
	dopplershift.from_prior(rng);

	pp = rng.rand();
	// this, too belongs to the noise process we're not using 
        noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
        noise_L = exp(log(1E-2*Data::get_instance().get_f_range())
                        + log(1E3)*rng.rand());

        calculate_mu();

}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;

	int which;

	// A 50-50 decision to perturb a 'multivariate' thing or a 'univariate' thing
	if(rng.rand() <= 0.5)
	{
		which = rng.rand_int(2);

		if(which == 0)	// Proposal for the RJObject describing the lines
			logH += dopplershift.perturb(rng);
		else			// Proposal for normals for the OU process
		{
			if(rng.rand() <= 0.5)	// Propose to move only one
			{
				which = rng.rand_int(noise_normals.size());
				logH -= -0.5*pow(noise_normals[which], 2);
				noise_normals[which] += rng.randh();
				logH += -0.5*pow(noise_normals[which], 2);
			}
			else					// AR(1) proposal
			{
				double theta = 2.*M_PI*pow(10., -6.*rng.rand());
				double cos_theta = cos(theta);
				double sin_theta = sin(theta);
				for(double& n: noise_normals)
					n = cos_theta*n + sin_theta*rng.randn();
			}
		}
	}
	else
	{
		which = rng.rand_int(4);
		if(which == 0)
		{
			background = log(background);
			background = (atan(background)/M_PI + 0.485)/0.97;
			background += rng.randh();
			background = mod(background, 1.);
			background = tan(M_PI*(0.97*background - 0.485));
			background = exp(background);
		}
		else if(which == 1)
		{
			pp += rng.randh();
			wrap(pp, 0., 1.);
		}
		else if(which == 2)
		{
			noise_sigma = log(noise_sigma);
			noise_sigma += log(1E3)*rng.randh();
			wrap(noise_sigma, log(1E-3), log(1.));
			noise_sigma = exp(noise_sigma);
		}
		else
		{
			noise_L = log(noise_L);
			noise_L += log(1E3)*rng.randh();
			wrap(noise_L, log(1E-2*Data::get_instance().get_f_range()), log(10.*Data::get_instance().get_f_range()));
			noise_L = exp(noise_L);
		}
	}

	calculate_mu();
	return logH;
}

double MyModel::log_likelihood() const
{
        const vector<double>& f_mid = data.get_f_mid();
        const vector<double>& y = data.get_y();
	const vector<double>& yerr = data.get_yerr();

        double logl = 0.;
	    for(size_t i=0; i<f_mid.size(); i++)
		{
			logl += -0.5*log(2.*M_PI) - log(yerr[i]) - pow(y[i]-mu[i], 2)/(2.*pow(yerr[i], 2));
 		}
	return logl;
}

void MyModel::print(std::ostream& out) const
{
        out<<background<<' ';
        dopplershift.print(out);

	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';

}

string MyModel::description() const
{
	return string("");
}

