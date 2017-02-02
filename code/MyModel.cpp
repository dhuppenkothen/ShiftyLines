#include "MyModel.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include "Lookup.h"
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace DNest4;

const Data& MyModel::data = Data::get_instance();
const int& nlines = Data::get_instance().get_nlines();
//const std::vector<double> line_pos = Data::get_instance().get_line_pos();
const PHAData& pha = Data::get_instance().get_pha();


MyModel::MyModel()
//:dopplershift(3*nlines+1, 5, false, MyConditionalPrior())
:dopplershift(3*nlines+1, 4, false, MyConditionalPrior())
,noise_normals(pha.bin_lo.size())
//,mu(data.get_pha().arf.specresp.size())
,mu(pha.counts.size())
{
}

double MyModel::gaussian_cdf(double x, double x0, double gamma)
{
	return 0.5*(1. + Lookup::erf((x-x0)/(gamma*sqrt(2.))));
}

template <typename ConstIntType, typename ConstFloatType,
            typename FloatType, typename IndexType, typename UIndexType>

void MyModel::rmf_fold(IndexType len_source, const ConstFloatType *source,
                IndexType len_num_groups, const ConstIntType *num_groups,
                IndexType len_first_chan, const ConstIntType *first_chan,
                IndexType len_num_chans, const ConstIntType *num_chans,
                IndexType len_response, const ConstFloatType *resp,
                IndexType len_counts, FloatType *counts,
                UIndexType offset)
{
    //int flag = 0;
//    if ( ( len_num_groups != len_source ) ||
//         ( len_first_chan != len_num_chans ))
//		throw RMFConvolutionFailure();

    // How many groups are in the current energy bin?
    IndexType current_num_groups = 0;

    // How many channels are in the current group?
    IndexType current_num_chans = 0;

    // What is the current channel of the output (counts) array?
    //register IndexType current_chan = 0;

    FloatType source_bin_ii;
    const ConstFloatType *resp_tmp = resp;
    const ConstIntType *first_chan_tmp = first_chan;
    const ConstIntType *num_chans_tmp = num_chans;
    FloatType *counts_tmp = counts;


    for (size_t ii = 0; ii < len_source; ii++ ) {

      // ii is the current energy bin
      source_bin_ii = source[ ii ];

      current_num_groups = num_groups[ ii ];

      while( current_num_groups ) {

//        if ( ( IndexType(first_chan_tmp - first_chan) >= len_num_chans ) ||
//             ( UIndexType(*first_chan_tmp) < offset ) )
//                //throw RMFConvolutionFailure();


        counts_tmp = counts + *first_chan_tmp - offset;
        current_num_chans = *num_chans_tmp;
        first_chan_tmp++;
        num_chans_tmp++;

//        if ( ( (IndexType(counts_tmp-counts) + current_num_chans) > len_counts )
//             ||
//             ( (IndexType(resp_tmp-resp) + current_num_chans) > len_response ) )
//                throw RMFConvolutionFailure();

        while ( current_num_chans ) {

          *counts_tmp += *resp_tmp * source_bin_ii;
          counts_tmp++;
          resp_tmp++;
          current_num_chans--;

        }
        current_num_groups--;

      }

    } // end for ii

  }



void MyModel::calculate_mu()
{ 
	// define line positions
	const vector<double>& line_pos = data.get_line_pos();

	// declare shifted line positions
	vector<double> line_pos_shifted;

	// get left and right boundaries of the wavelength bins
        //const vector<double>& f_left = data.get_f_left();
        //const vector<double>& f_right = data.get_f_right();

        // NEW VERSION: get PHA data from FITS file:
        // doesn't do anything yet, just testing whether I can load the data!
	// NOTE: Just having this line in the code (uncommented, of course), makes
	// the code slower by at least a factor of 3! Not sure why that is, but I 
	// should probably figure that out
        const PHAData& pha = data.get_pha();

	const vector<double>& f_left = pha.bin_lo;
	const vector<double>& f_right = pha.bin_hi;

	// read out the ARF
        //const vector<double>& _arf = pha.arf.specresp;
	//const RMFData& rmf = pha.rmf;

	// assign constant background to model
	mu.assign(mu.size(), background); //old version
 
	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();
 
//	vector<double> amplitude, logq, sign, width;
	vector<double> amplitude, sign, width;
	double dshift;

	vector<double> mu_temp(mu.size(), 0.0);

        for(size_t j=0; j<dopplershiftcomponents.size(); j++)
        {
	        // assign 0 to all shifted line positions
	        line_pos_shifted.assign(line_pos.size(), 0.0);

		amplitude.assign(nlines, 0.);
//		logq.assign(nlines, 0.);
		sign.assign(nlines, 0.);
		width.assign(nlines, 0.);

		dshift = dopplershiftcomponents[j][0];
	
		for (int i=0; i<nlines; i++)
			{
				line_pos_shifted[i] = line_pos[i]*(1. + dshift);
				amplitude[i] = exp(dopplershiftcomponents[j][i+1]);
				//logq[i] = dopplershiftcomponents[j][i+1+nlines];
				sign[i] = dopplershiftcomponents[j][i+1+2*nlines];		
				width[i] = exp(dopplershiftcomponents[j][i+1+nlines]);
			}
	
		
		int s=0;	
                for(size_t i=0; i<mu.size(); i++)
                {
			for (int k=0; k<nlines; k++)
				{
			// Integral over the Lorentzian distribution
					if (sign[k] < dopplershift.get_conditional_prior().get_pp()) 
						s = -1;
					else 
						s = 1;
 
				//	mu_temp[i] += s*amplitude[k]/(width[k]*sqrt(2.*M_PI))*exp(-pow(f_mid[i]-line_pos_shifted[k],2)/(2.*pow(width[k],2)));
					mu_temp[i] += s*amplitude[k]*(gaussian_cdf(f_right[i], line_pos_shifted[k], width[k])
								- gaussian_cdf(f_left[i], line_pos_shifted[k], width[k]));

				}

                } 

	}


	/////// NOTE ////////////////////
	// the code below cuts out a small part of the original spectrum for testing purposes. 
	// NEED TO DELETE THIS ONCE EVERYTHING WORKS!!!!!
        //

        const vector<double>& f_mid_old = data.get_f_mid();
        const vector<double>& f_mid = pha.bin_lo;


        // I'm only interested in a specific region of the spectrum
        // right now, so let's only look at that!
        double f_min = f_mid_old[0];
        double f_max = f_mid_old[f_mid_old.size()-1];


        for (size_t k=0; k<pha.bin_lo.size(); k++)
                {
                        if (f_mid[k] < f_min)
                                continue;
                        if (f_mid[k] > f_max)
                                continue;
                        else
                                {
                                        mu_small.push_back(mu[k]);
                                }

                }



        // fold through the ARF
        // code taken from sherpa
        for (size_t ii = 0; ii < mu.size(); ii++ )
             mu[ ii ] *= pha.arf.specresp[ ii ];


        // Compute the OU process
        vector<double> y(mu.size());
        double alpha = exp(-1./noise_L);

	// noise process could come both from the source or the detector!
	// which is why I put it in between the ARF and the RMF
        for(size_t i=0; i<mu.size(); i++)
        {
                if(i==0)
                        y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
                else
                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
                mu[i] *= exp(y[i]);
        }


	// Remove all counts < zero, because that would be just silly!
        
	for (size_t i=0; i<mu.size(); i++)
		{
//			mu[i]  += mu_temp[i];//exp(mu_temp[i]);
			if (mu[i] < 0.0)
				mu[i] = 0.0;
		}

        counts.assign(mu.size(), 0.0);


        rmf_fold(mu.size(), &mu[0], 
		 pha.rmf.n_grp.size(), &pha.rmf.n_grp[0],
             	 pha.rmf.f_chan.size(), &pha.rmf.f_chan[0],
		 pha.rmf.n_chan.size(), &pha.rmf.n_chan[0],
		 pha.rmf.matrix.size(), &pha.rmf.matrix[0],
		 counts.size(), &counts[0],
		 pha.rmf.offset);


}

void MyModel::from_prior(RNG& rng)
{
	background = tan(M_PI*(0.97*rng.rand() - 0.485));
	background = exp(background);
	dopplershift.from_prior(rng);

	// this, too belongs to the noise process we're not using 
        noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
        noise_L = exp(log(0.1*Data::get_instance().get_f_range())
                        + log(1E2)*rng.rand());

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
		which = rng.rand_int(3);
		if(which == 0)
		{
			background = log(background);
			background = (atan(background)/M_PI + 0.485)/0.97;
			background += rng.randh();
			background = mod(background, 1.);
			background = tan(M_PI*(0.97*background - 0.485));
			background = exp(background);
		}
		if(which == 1)
		{
			noise_sigma = log(noise_sigma);
			noise_sigma += log(1E3)*rng.randh();
			wrap(noise_sigma, log(1E-3), log(1.));
			noise_sigma = exp(noise_sigma);
		}
		else
		{
			noise_L = log(noise_L);
			noise_L += log(1E2)*rng.randh();
			wrap(noise_L, log(0.1*Data::get_instance().get_f_range()), log(10.*Data::get_instance().get_f_range()));
			noise_L = exp(noise_L);
		}
	}

	calculate_mu();
	return logH;
}

double MyModel::log_likelihood() const
{
        const vector<double>& f_mid_old = data.get_f_mid();

//        const vector<double>& y = data.get_y();
//	const vector<double>& yerr = data.get_yerr();
	const PHAData& pha = data.get_pha();
	const vector<double>& y = pha.counts;
	const vector<double>& f_mid = pha.bin_lo;

	// Very dumb approximation for error
	// NEED TO IMPLEMENT POISSON LIKELIHOOD!!!
        vector<double> yerr(y.size(), 0.0);

	for (size_t j=0; j<y.size(); j++)
		{
			yerr[j] += sqrt(y[j]);
		} 



	// I'm only interested in a specific region of the spectrum
	// right now, so let's only look at that!
	double f_min = f_mid_old[0];
	double f_max = f_mid_old[f_mid_old.size()-1];

	vector<double> y_small, yerr_small, mu_small;
	for (size_t k=0; k<y.size(); k++)
		{
			if (f_mid[k] < f_min)
				continue;
			if (f_mid[k] > f_max)
				continue;
			else
				{
					y_small.push_back(y[k]);
					yerr_small.push_back(yerr[k]);
					mu_small.push_back(mu[k]);
				} 

		}

        double logl = 0.;
	    for(size_t i=0; i<y_small.size(); i++)
		{
			logl += -mu_small[i] + y_small[i]*log(mu_small[i]) - gsl_sf_lngamma(y_small[i] + 1.);

			//logl += -0.5*log(2.*M_PI) - log(yerr_small[i]) - pow(y_small[i]-mu_small[i], 2)/(2.*pow(yerr_small[i], 2));
 		}
	return logl;
}

void MyModel::print(std::ostream& out) const
{
        out<<background<<' '<<noise_L<<' '<<noise_sigma<<' ';
        dopplershift.print(out);

	for(size_t i=0; i<mu.size(); i++)
		out<<mu_small[i]<<' ';

}

string MyModel::description() const
{
	return string("");
}

