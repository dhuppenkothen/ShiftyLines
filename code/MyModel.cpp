#include "MyModel.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest4;

const Data& MyModel::data = Data::get_instance();
#include <iostream>

const int nlines = Data::get_instance().get_nlines();

MyModel::MyModel()
:dopplershift(3*nlines+1, 5, false, MyConditionalPrior())
//:dopplershift(25, 5, false, MyConditionalPrior())
,mu(data.get_f_left().size())
{
}

double MyModel::gaussian_cdf(double x, double x0, double gamma)
{
	return 0.5*(1. + erf((x-x0)/(gamma*pow(2., 0.5))));
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

	// assign constant background to model
	mu.assign(mu.size(), exp(background));

//	mu.assign(mu.size(), background); //old version
 
	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();

	vector<double> amplitude, logq, sign, width;
	double dshift;

	vector<double> mu_temp;
	mu_temp.assign(mu.size(), 0.);

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
				line_pos_shifted[i] = line_pos[0] + dshift;
				amplitude[i] = exp(dopplershiftcomponents[j][i+1]);
				logq[i] = dopplershiftcomponents[j][i+1+nlines];
				sign[i] = dopplershiftcomponents[j][i+1+2*nlines]	;		
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
 
					mu_temp[i] += s*amplitude[k]*(gaussian_cdf(f_right[i], line_pos_shifted[k], width[k])
								- gaussian_cdf(f_left[i], line_pos_shifted[k], width[k]));

				}

                } 

	// this is a OU process; we're currently not going to use that, but it might come 
	// in handy later, so we'll leave it commented out at the appropriate places
//        vector<double> y(mu.size());
//        double alpha = exp(-1./noise_L);
 
//        for(size_t i=0; i<mu.size(); i++)
//        {
//                if(i==0)
//                        y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
//                else
//                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
//                mu[i] *= exp(y[i]);
//        }


        }
	for (size_t i=0; i<mu.size(); i++)
		mu[i]  += exp(mu_temp[i]);
}

void MyModel::from_prior(RNG& rng)
{
	background = tan(M_PI*(0.97*rng.rand() - 0.485));
	background = exp(background);
	dopplershift.from_prior(rng);

	pp = rng.rand();
	// this, too belongs to the noise process we're not using 
//        noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
//        noise_L = exp(log(1E-2*Data::get_instance().get_t_range())
//                        + log(1E3)*rng.rand());

        calculate_mu();

}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;

        if(rng.rand() <= 0.2)
        {
//                for(size_t i=0; i<mu.size(); i++)
//                        mu[i] -= background;
                background = log(background);
                background = (atan(background)/M_PI + 0.485)/0.97;
                background += pow(10., 1.5 - 6.*rng.rand())*rng.randn();
                background = mod(background, 1.);
                background = tan(M_PI*(0.97*background - 0.485));
                background = exp(background);

			calculate_mu();

//                for(size_t i=0; i<mu.size(); i++)
//                        mu[i] += background;
        }
	else if(rng.rand() <= 0.2)
	{
		pp += rng.randh();
		wrap(pp, 0., 1.);
	}
	else
	{
		logH += dopplershift.perturb(rng);
		calculate_mu();
	}
//        else if(rng.rand() <= 0.5)
//        {
//                noise_sigma = log(noise_sigma);
//                noise_sigma += log(1E3)*rng.randh();
//                wrap(noise_sigma, log(1E-3), log(1.));
//                noise_sigma = exp(noise_sigma);

//                noise_L = log(noise_L);
//                noise_L += log(1E3)*rng.randh();
//                wrap(noise_L, log(1E-2*Data::get_instance().get_t_range()), log(10.*Data::get_instance().get_t_range()));
//                noise_L = exp(noise_L);

//                calculate_mu();
//        }
//        else
//        {
//                int num = exp(log((double)noise_normals.size())*rng.rand());
//                for(int i=0; i<num; i++)
//                {
//                        int k = rng.rand_int(noise_normals.size());
//                        noise_normals[k] = rng.randn();
//                }
//                calculate_mu();
//        }


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
			logl += -0.5*log(2.*M_PI) - 0.5*log(yerr[i]) - pow(y[i]-mu[i], 2)/(2.*pow(yerr[i], 2));
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

