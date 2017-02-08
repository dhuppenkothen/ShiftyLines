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
const PHAData& pha_heg_p1 = Data::get_instance().get_pha_heg_p1();
const PHAData& pha_heg_m1 = Data::get_instance().get_pha_heg_m1();
//const PHAData& pha_meg_p1 = Data::get_instance().get_pha_meg_p1();
//const PHAData& pha_meg_m1 = Data::get_instance().get_pha_meg_m1();


MyModel::MyModel()
//:dopplershift(3*nlines+1, 5, false, MyConditionalPrior())
:dopplershift(3*nlines+1, 4, false, MyConditionalPrior())
,noise_normals(pha_heg_p1.bin_lo.size())
//,mu(data.get_pha().arf.specresp.size())
,mu1(pha_heg_p1.counts.size())
,mu2(pha_heg_m1.counts.size())
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
	vector<double> line_pos_shifted(line_pos.size());

	// get left and right boundaries of the wavelength bins
        //const vector<double>& f_left = data.get_f_left();
        //const vector<double>& f_right = data.get_f_right();

        // NEW VERSION: get PHA data from FITS file:
        // doesn't do anything yet, just testing whether I can load the data!
	// NOTE: Just having this line in the code (uncommented, of course), makes
	// the code slower by at least a factor of 3! Not sure why that is, but I 
	// should probably figure that out

	const vector<double>& f_left = pha_heg_p1.bin_lo;
	const vector<double>& f_right = pha_heg_p1.bin_hi;

	// assign constant background to model
	mu1.assign(mu1.size(), background);
        mu2.assign(mu2.size(), background); 

	mu.assign(mu1.size(), 0.0); // array 

	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();
 
	vector<double> amplitude(nlines), sign(nlines), width(nlines);
	double dshift;

//        // I'm only interested in a specific region of the spectrum
//        // right now, so let's only look at that!
	const double& f_min = data.get_f_min();
        const double& f_max = data.get_f_max();

        for(size_t j=0; j<dopplershiftcomponents.size(); j++)
        {

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
        	        for(size_t i=0; i<mu1.size(); i++)
        	        {
		                if (f_left[i] < f_min)                        
		                       continue;
		                if (f_right[i] > f_max)
	        	               continue;
	       		        else
	                	{       


					for (int k=0; k<nlines; k++)
						{
						// Integral over the Lorentzian distribution
						if (sign[k] < dopplershift.get_conditional_prior().get_pp()) 
							s = -1;
						else 
							s = 1;
						if ((std::abs(f_right[i] - line_pos_shifted[k]) < 5.*width[k]) && 
						   (std::abs(f_left[i] - line_pos_shifted[k]) < 5.*width[k])) 
							mu[i] += s*amplitude[k]*(gaussian_cdf(f_right[i], line_pos_shifted[k], width[k])
										- gaussian_cdf(f_left[i], line_pos_shifted[k], width[k]));
						}
                		}	 
			}
	}

   
        // fold through the ARF
        // code taken from sherpa
        for (size_t ii = 0; ii < mu1.size(); ii++ )
		{
			mu1[ ii ] = (mu1[ ii ] + mu[ ii ]) * pha_heg_p1.arf.specresp[ ii ];
             		mu2[ ii ] = inst_fac * (mu2[ ii ] + mu[ ii]) * pha_heg_m1.arf.specresp[ ii ];
		}

        // Compute the OU process
        vector<double> y(mu1.size());
        double alpha = exp(-1./noise_L);

	// noise process could come both from the source or the detector!
 	// which is why I put it in between the ARF and the RMF
	for(size_t i=0; i<mu1.size(); i++)
	{
		if (f_left[i] < f_min)
			y[i]=1.0;
		if (f_right[i] > f_max)
             		y[i]=1.0;

	        if((f_left[i] < f_min) && (f_right[i] > f_min ))
	                y[i] = noise_sigma/sqrt(1. - alpha*alpha)*noise_normals[i];
	        else
	                y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];
	        mu1[i] *= exp(y[i]);
                mu2[i] *= exp(y[i]);

	}

        counts1.assign(mu1.size(), 0.0);
        counts2.assign(mu2.size(), 0.0);


        rmf_fold(mu1.size(), &mu1[0], 
	  pha_heg_p1.rmf.n_grp.size(), &pha_heg_p1.rmf.n_grp[0],
	  pha_heg_p1.rmf.f_chan.size(), &pha_heg_p1.rmf.f_chan[0],
	  pha_heg_p1.rmf.n_chan.size(), &pha_heg_p1.rmf.n_chan[0],
	  pha_heg_p1.rmf.matrix.size(), &pha_heg_p1.rmf.matrix[0],
	  counts1.size(), &counts1[0],
	  pha_heg_p1.rmf.offset);

        rmf_fold(mu2.size(), &mu2[0],
          pha_heg_m1.rmf.n_grp.size(), &pha_heg_m1.rmf.n_grp[0],
          pha_heg_m1.rmf.f_chan.size(), &pha_heg_m1.rmf.f_chan[0],
          pha_heg_m1.rmf.n_chan.size(), &pha_heg_m1.rmf.n_chan[0],
          pha_heg_m1.rmf.matrix.size(), &pha_heg_m1.rmf.matrix[0],
          counts2.size(), &counts2[0],
          pha_heg_m1.rmf.offset);


	mu1.resize(data.get_pha_heg_p1().bin_lo.size());
        mu2.resize(data.get_pha_heg_m1().bin_lo.size());

}

void MyModel::from_prior(RNG& rng)
{
	background = tan(M_PI*(0.97*rng.rand() - 0.485));
	background = exp(background);
        inst_fac = tan(M_PI*(0.97*rng.rand() - 0.485));
        inst_fac = exp(inst_fac);

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
                if(which == 1)
                {
                        inst_fac = log(inst_fac);
                        inst_fac = (atan(inst_fac)/M_PI + 0.485)/0.97;
                        inst_fac += rng.randh();
                        inst_fac = mod(inst_fac, 1.);
                        inst_fac = tan(M_PI*(0.97*inst_fac - 0.485));
                        inst_fac = exp(inst_fac);
                }

		if(which == 2)
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
	const PHAData& pha_heg_p1 = data.get_pha_heg_p1();
        const PHAData& pha_heg_m1 = data.get_pha_heg_m1();

	const vector<double>& y1 = pha_heg_p1.counts;
	const vector<double>& y2 = pha_heg_m1.counts;
	const vector<double>& f_left = pha_heg_p1.bin_lo;
        const vector<double>& f_right = pha_heg_p1.bin_hi;

	// I'm only interested in a specific region of the spectrum
	// right now, so let's only look at that!

        const double& f_min = data.get_f_min();
        const double& f_max = data.get_f_max();

        double logl1 = 0.;
	double logl2 = 0.;
	    for(size_t i=0; i<y1.size(); i++)
		{
                        if (f_left[i] < f_min)
                                continue;
                        if (f_right[i] > f_max)
                                continue;
			else
				{
					logl1 += -mu1[i] + y1[i]*log(mu1[i]) - gsl_sf_lngamma(y1[i] + 1.);
                                        logl2 += -mu2[i] + y2[i]*log(mu2[i]) - gsl_sf_lngamma(y2[i] + 1.);
				}
 		}
	return logl1+logl2;
}

void MyModel::print(std::ostream& out) const
{
        out<<background<<' '<<inst_fac<<' '<<noise_L<<' '<<noise_sigma<<' ';
        dopplershift.print(out);

	for(size_t i=0; i<mu1.size(); i++)
		out<<mu1[i]<<' ';

        for(size_t i=0; i<mu2.size(); i++)
                out<<mu2[i]<<' ';

  
}

string MyModel::description() const
{
	return string("");
}

