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
const PHAData& pha_heg_p = Data::get_instance().get_pha_heg_p();
//const PHAData& pha_heg_m = Data::get_instance().get_pha_heg_m();
//const PHAData& pha_meg_p = Data::get_instance().get_pha_meg_p();
//const PHAData& pha_meg_m = Data::get_instance().get_pha_meg_m();

// Initialise the static distribution
const DNest4::Cauchy MyModel::cauchy(0.0, 1.0);

MyModel::MyModel()
:dopplershift(3*nlines+1, 4, false, MyConditionalPrior())
,noise_normals_h(pha_heg_p.bin_lo.size())
//,noise_normals_m(pha_meg_p.bin_lo.size())
,mu_hp(pha_heg_p.counts.size())
//,mu_hm(pha_heg_m.counts.size())
//,mu_mp(pha_meg_p.counts.size())
//,mu_mm(pha_meg_m.counts.size())
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

        // NEW VERSION: get PHA data from FITS file:
        // doesn't do anything yet, just testing whether I can load the data!
	// NOTE: Just having this line in the code (uncommented, of course), makes
	// the code slower by at least a factor of 3! Not sure why that is, but I 
	// should probably figure that out

	const vector<double>& f_left_h = pha_heg_p.bin_lo;
	const vector<double>& f_right_h = pha_heg_p.bin_hi;
//        const vector<double>& f_left_m = pha_meg_p.bin_lo;
//        const vector<double>& f_right_m = pha_meg_p.bin_hi;

	const string& eunit = pha_heg_p.bin_lo_unit;

	// assign constant background to mode
	mu_hp.assign(mu_hp.size(), 0.0);
//        mu_hm.assign(mu_hm.size(), 0.0); 
//        mu_mp.assign(mu_mp.size(), 0.0);
//        mu_mm.assign(mu_mm.size(), 0.0);

	mu_h.assign(mu_hp.size(), 0.0); // array 
//        mu_m.assign(mu_mp.size(), 0.0); // array 

        mu_hp_out.assign(mu_hp.size(), 0.0);
	mu_hp_specresp.assign(mu_hp.size(), 0.0);
//        mu_hm_out.assign(mu_hm.size(), 0.0);
//        mu_mp_out.assign(mu_mp.size(), 0.0);
//        mu_mm_out.assign(mu_mm.size(), 0.0);


	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();
 
	vector<double> amplitude(nlines), sign(nlines), width(nlines);
	double dshift;

//        // I'm only interested in a specific region of the spectrum
//        // right now, so let's only look at that!
	const double& f_min = data.get_f_min();
	const double f_max = data.get_f_max();

        for(size_t j=0; j<dopplershiftcomponents.size(); j++)
        {

			dshift = dopplershiftcomponents[j][0];
	
			for (int i=0; i<nlines; i++)
				{
					if (eunit == "keV")
						{
						line_pos_shifted[i] = line_pos[i]*(1.0/(1. + dshift));
						}
					if (eunit == "angstrom")
						line_pos_shifted[i] = line_pos[i]*(1. + dshift);
					amplitude[i] = exp(dopplershiftcomponents[j][i+1]);
					//logq[i] = dopplershiftcomponents[j][i+1+nlines];
      					sign[i] = dopplershiftcomponents[j][i+1+2*nlines];		
					width[i] = exp(dopplershiftcomponents[j][i+1+nlines]);
				}
		
			int sh=0;	
        	        for(size_t i=0; i<mu_h.size(); i++)
        	        {
		                if (f_left_h[i] < f_min)                        
		                       continue;
		                if (f_right_h[i] > f_max)
	        	               continue;
	       		        else
	                	{       


					for (int k=0; k<nlines; k++)
						{
						// Integral over the Lorentzian distribution
						if (sign[k] < dopplershift.get_conditional_prior().get_pp()) 
							sh = -1;
						else 
							sh = 1;
						if ((std::abs(f_right_h[i] - line_pos_shifted[k]) < 5.*width[k]) && 
						   (std::abs(f_left_h[i] - line_pos_shifted[k]) < 5.*width[k])) 
							mu_h[i] += sh*amplitude[k]*(gaussian_cdf(f_right_h[i], line_pos_shifted[k], width[k])
										- gaussian_cdf(f_left_h[i], line_pos_shifted[k], width[k]));
					}
                		}	 
			}

//                        int sm=0;
//                        for(size_t i=0; i<mu_m.size(); i++)
//                        {
//                                if (f_left_m[i] < f_min)
//                                       continue;
//                                if (f_right_m[i] > f_max)
//                                       continue;
//                                else
//                                {
//
//
//                                        for (int k=0; k<nlines; k++)
//                                                {
//                                                // Integral over the Lorentzian distribution
//                                                if (sign[k] < dopplershift.get_conditional_prior().get_pp())
//                                                        sm = -1;
//                                                else
//                                                        sm = 1;
//                                                if ((std::abs(f_right_m[i] - line_pos_shifted[k]) < 5.*width[k]) &&
//                                                   (std::abs(f_left_m[i] - line_pos_shifted[k]) < 5.*width[k]))
//                                                        mu_m[i] += sm*amplitude[k]*(gaussian_cdf(f_right_m[i], line_pos_shifted[k], width[k])
//                                                                                - gaussian_cdf(f_left_m[i], line_pos_shifted[k], width[k]));
//                                                }
//                                }
//                        }
//
	}
//
//
   
        // fold through the ARF
        // code taken from sherpa
        for (size_t ii = 0; ii < mu_h.size(); ii++ )
		{
			mu_hp_out[ ii ] =  exp(log(background) + mu_h[ ii ]);
//			mu_hm_out[ ii ] =  (mu_hm[ ii ] + mu_h[ ii]);

			mu_hp[ ii ] = exp(log(background) + mu_h[ ii ]);
                        mu_hp[ ii ] *= pha_heg_p.arf.specresp[ ii ];

			mu_hp_specresp[ ii ] = (mu_hp[ ii ] + mu_h[ ii ]) * pha_heg_p.arf.specresp[ ii ]; 
//            		mu_hm[ ii ] = mu_hm_out[ ii ] * pha_heg_m.arf.specresp[ ii ];

		}
 
//        for (size_t ii = 0; ii < mu_m.size(); ii++ )
//                {
//			mu_mp_out[ ii ] = exp(log(background) + mu_m[ ii ]);
//			mu_mm_out[ ii ] = exp(log(background) + mu_m[ ii]);
//
//                        mu_mp[ ii ] = mu_mp_out[ ii ] * pha_meg_p.arf.specresp[ ii ];
//                        mu_mm[ ii ] = mu_mm_out[ ii ] * pha_meg_m.arf.specresp[ ii ];
//
//                }
//


        // Compute the OU process
        //vector<double> y_h(mu_h.size());
        double alpha = exp(-1./noise_L);

        //vector<double> y_m(mu_m.size());

        y_h.assign(mu_hp.size(), 0.0);
//        y_m.assign(mu_hm.size(), 0.0);

	// noise process could come both from the source or the detector!
 	// which is why I put it in between the ARF and the RMF
	for(size_t i=0; i<mu_h.size(); i++)
	{
	        if(i == 0)
        	    	y_h[i] = noise_sigma*noise_normals_h[i];
       	 	else
            		y_h[i] = alpha*y_h[i-1] + noise_sigma*noise_normals_h[i];

	        mu_hp[i] *= exp(y_h[i]);
		//mu_hm[i] *= (inst_fac_hm * exp(y_h[i]));
		//mu_hp_out[ i ] *= exp(y_h[ i ]);
                //mu_hm_out[ i ] *= exp(y_h[ i ]);

	}

        counts_hp.assign(mu_hp.size(), 0.0);
//        counts_hm.assign(mu_hm.size(), 0.0);

        // noise process could come both from the source or the detector!
        // which is why I put it in between the ARF and the RMF
//        for(size_t i=0; i<mu_m.size(); i++)
//        {
//                if(i == 0)
//                        y_m[i] = noise_sigma*noise_normals_m[i];
//                else
//                        y_m[i] = alpha*y_m[i-1] + noise_sigma*noise_normals_m[i];

//                mu_mp[i] *= (inst_fac_mp * exp(y_m[i]));
 //               mu_mm[i] *= (inst_fac_mm * exp(y_m[i]));
 //               mu_mp_out[ i ] *= exp(y_m[ i ]);
 //               mu_mm_out[ i ] *= exp(y_m[ i ]);

//        }

//        counts_mp.assign(mu_mp.size(), 0.0);
//        counts_mm.assign(mu_mm.size(), 0.0);

        rmf_fold(mu_hp.size(), &mu_hp[0], 
	  pha_heg_p.rmf.n_grp.size(), &pha_heg_p.rmf.n_grp[0],
	  pha_heg_p.rmf.f_chan.size(), &pha_heg_p.rmf.f_chan[0],
	  pha_heg_p.rmf.n_chan.size(), &pha_heg_p.rmf.n_chan[0],
	  pha_heg_p.rmf.matrix.size(), &pha_heg_p.rmf.matrix[0],
	  counts_hp.size(), &counts_hp[0],
	  pha_heg_p.rmf.offset);

//        rmf_fold(mu_hm.size(), &mu_hm[0],
//          pha_heg_m.rmf.n_grp.size(), &pha_heg_m.rmf.n_grp[0],
//          pha_heg_m.rmf.f_chan.size(), &pha_heg_m.rmf.f_chan[0],
//          pha_heg_m.rmf.n_chan.size(), &pha_heg_m.rmf.n_chan[0],
//          pha_heg_m.rmf.matrix.size(), &pha_heg_m.rmf.matrix[0],
//          counts_hm.size(), &counts_hm[0],
//          pha_heg_m.rmf.offset);
//
//        rmf_fold(mu_mp.size(), &mu_mp[0],
//          pha_meg_p.rmf.n_grp.size(), &pha_meg_p.rmf.n_grp[0],
//          pha_meg_p.rmf.f_chan.size(), &pha_meg_p.rmf.f_chan[0],
//          pha_meg_p.rmf.n_chan.size(), &pha_meg_p.rmf.n_chan[0],
//          pha_meg_p.rmf.matrix.size(), &pha_meg_p.rmf.matrix[0],
//          counts_mp.size(), &counts_mp[0],
//          pha_meg_p.rmf.offset);
//
//        rmf_fold(mu_mm.size(), &mu_mm[0],
//          pha_meg_m.rmf.n_grp.size(), &pha_meg_m.rmf.n_grp[0],
//          pha_meg_m.rmf.f_chan.size(), &pha_meg_m.rmf.f_chan[0],
//          pha_meg_m.rmf.n_chan.size(), &pha_meg_m.rmf.n_chan[0],
//          pha_meg_m.rmf.matrix.size(), &pha_meg_m.rmf.matrix[0],
//          counts_mm.size(), &counts_mm[0],
//          pha_meg_m.rmf.offset);
//
        counts_hp.resize(data.get_pha_heg_p().bin_lo.size());
//        counts_hm.resize(data.get_pha_heg_m().bin_lo.size());
//
//	counts_mp.resize(data.get_pha_meg_p().bin_lo.size());
//        counts_mm.resize(data.get_pha_meg_m().bin_lo.size());
//
}

void MyModel::from_prior(RNG& rng)
{
    do
    {
    	background = cauchy.generate(rng);
    }while(std::abs(background) > 25.0);
    background = exp(background);


//	do
//	{
//		inst_fac_hm = cauchy.generate(rng);
//	}while(std::abs(inst_fac_hm) > 25.0);
//    	inst_fac_hm = exp(inst_fac_hm);
//
//
//        do
//        {
//                inst_fac_mp = cauchy.generate(rng);
//        }while(std::abs(inst_fac_mp) > 25.0);
//        inst_fac_mp = exp(inst_fac_mp);
//
//        do
//        {
//                inst_fac_mm = cauchy.generate(rng);
//        }while(std::abs(inst_fac_mm) > 25.0);
//        inst_fac_mm = exp(inst_fac_mm);

	dopplershift.from_prior(rng);

	// this, too belongs to the noise process we're not using 
        noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
        noise_L = exp(log(0.01*Data::get_instance().get_f_range())
                        + log(1000)*rng.rand());

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
				int i = rng.rand_int(noise_normals_h.size());
				logH -= -0.5*pow(noise_normals_h[i], 2);
				noise_normals_h[i] += rng.randh();
				logH += -0.5*pow(noise_normals_h[i], 2);

 //                               i = rng.rand_int(noise_normals_m.size());
 //                               logH -= -0.5*pow(noise_normals_m[i], 2);
 //                               noise_normals_m[i] += rng.randh();
 //                               logH += -0.5*pow(noise_normals_m[i], 2);



			}
			else					// Regenerate many
			{
                int reps = (int)pow(noise_normals_h.size(), rng.rand());
                for(int i=0; i<reps; ++i)
                {
                    int k = rng.rand_int(noise_normals_h.size());
                    noise_normals_h[k] = rng.randn();
                }

//                reps = (int)pow(noise_normals_m.size(), rng.rand());
//                for(int i=0; i<reps; ++i)
//                {
//                    int k = rng.rand_int(noise_normals_m.size());
//                    noise_normals_m[k] = rng.randn();
//                }
			}
		}
	}
	else
	{
		which = rng.rand_int(3);
		if(which == 0)
		{
            		background = log(background);
            		logH += cauchy.perturb(background, rng);
            		if(std::abs(background) > 25.0)
                		logH = -1E300;
            		background = exp(background);
		}
//                else if(which == 1)
//                {
//                        inst_fac_hm = log(inst_fac_hm);
//                        logH += cauchy.perturb(inst_fac_hm, rng);
//                        if(std::abs(inst_fac_hm) > 25.0)
//                                logH = -1E300;
//			inst_fac_hm = exp(inst_fac_hm);
//                }
//
//                else if(which == 2)
//                {
//			inst_fac_mp = log(inst_fac_mp);
//                        logH += cauchy.perturb(inst_fac_mp, rng);
//                        if(std::abs(inst_fac_mp) > 25.0)
//                                logH = -1E300;
//			inst_fac_mp = exp(inst_fac_mp);
//                }
//
//                else if(which == 3)
//                {
//			inst_fac_mm = log(inst_fac_mm);
//                        logH += cauchy.perturb(inst_fac_mm, rng);
//                        if(std::abs(inst_fac_mm) > 25.0)
//                                logH = -1E300;
//			inst_fac_mm = exp(inst_fac_mm);
//                }
//

		else if(which == 1)
		{
                        noise_sigma = log(noise_sigma);
                        logH += cauchy.perturb(noise_sigma, rng);
                        if(noise_sigma < -20.0 || noise_sigma > 0.0)
                                return -1E300;
                        noise_sigma = exp(noise_sigma);

		}
		else
		{
			noise_L = log(noise_L);
			noise_L += log(1E5)*rng.randh();
			wrap(noise_L, log(0.01*Data::get_instance().get_f_range()), log(1000.*Data::get_instance().get_f_range()));
			noise_L = exp(noise_L);
		}
	}

	calculate_mu();
	return logH;
}

double MyModel::log_likelihood() const
{


	const PHAData& pha_heg_p = data.get_pha_heg_p();
//        const PHAData& pha_heg_m = data.get_pha_heg_m();
//        const PHAData& pha_meg_p = data.get_pha_meg_p();
//        const PHAData& pha_meg_m = data.get_pha_meg_m();


	const vector<double>& y_hp = pha_heg_p.counts;
//	const vector<double>& y_hm = pha_heg_m.counts;
	const vector<double>& f_left_h = pha_heg_p.bin_lo;
        const vector<double>& f_right_h = pha_heg_p.bin_hi;

//        const vector<double>& y_mp = pha_meg_p.counts;
//        const vector<double>& y_mm = pha_meg_m.counts;
//        const vector<double>& f_left_m = pha_meg_p.bin_lo;
//        const vector<double>& f_right_m = pha_meg_p.bin_hi;

	// I'm only interested in a specific region of the spectrum
	// right now, so let's only look at that!

        const double& f_min = data.get_f_min();
        //const double& f_max = data.get_f_max();
	const double& f_max = data.get_f_max();

        double logl_hp = 0.;
//	double logl_hm = 0.;
	    for(size_t i=0; i<y_hp.size(); i++)
		{
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
			else
				{
					logl_hp += -counts_hp[i] + y_hp[i]*log(counts_hp[i]) - gsl_sf_lngamma(y_hp[i] + 1.);
//                                        logl_hm += -counts_hm[i] + y_hm[i]*log(counts_hm[i]) - gsl_sf_lngamma(y_hm[i] + 1.);

				}
 		}

//        double logl_mp = 0.;
//        double logl_mm = 0.;
//            for(size_t i=0; i<y_mp.size(); i++)
//                {
//                        if (f_left_m[i] < f_min)
//                                continue;
//                        if (f_right_m[i] > f_max)
//                                continue;
//                        else
//                                {
//                                        logl_mp += -counts_mp[i] + y_mp[i]*log(counts_mp[i]) - gsl_sf_lngamma(y_mp[i] + 1.);
//                                        logl_mm += -counts_mm[i] + y_mm[i]*log(counts_mm[i]) - gsl_sf_lngamma(y_mm[i] + 1.);
//                                }
//                }
//


	return logl_hp;
//	return logl_hp + logl_hm + logl_mp + logl_mm;
}

void MyModel::print(std::ostream& out) const
{
//        out<<background<<' '<<inst_fac_hm<<' '<<inst_fac_mp<<' '<<inst_fac_mm<<' '<<noise_L<<' '<<noise_sigma<<' ';
	out<<background<<' '<<noise_L<<' '<<noise_sigma<<' ';

        dopplershift.print(out);

        const double& f_min = data.get_f_min();
        const double& f_max = data.get_f_max();
        const vector<double>& f_left_h = pha_heg_p.bin_lo;
        const vector<double>& f_right_h = pha_heg_p.bin_hi;

//        const vector<double>& f_left_m = pha_meg_p.bin_lo;
//        const vector<double>& f_right_m = pha_meg_p.bin_hi;

        for(size_t i=0; i<mu_hp_out.size(); i++)
                {
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
                        else
                                out<<y_h[i]<<' ';
		}


//        for(size_t i=0; i<mu_hp_out.size(); i++)
//                {
//                        if (f_left_h[i] < f_min)
//                                continue;
//                        if (f_right_h[i] > f_max)
//                                continue;
//                        else
//                                out<<y_m[i]<<' ';
//		}

        for(size_t i=0; i<mu_hp_out.size(); i++)
                {
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
                        else
                                out<<mu_h[i]<<' ';
                }


        for(size_t i=0; i<mu_hp_out.size(); i++)
                {
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
                        else
                                out<<mu_hp_out[i]<<' ';
                }
	

        for(size_t i=0; i<mu_hp_out.size(); i++)
                {
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
                        else
                                out<<mu_hp_specresp[i]<<' ';
                }


        for(size_t i=0; i<mu_hp_out.size(); i++)
                {
                        if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
                        else
                                out<<mu_hp[i]<<' ';
                }


//        for(size_t i=0; i<mu_hm_out.size(); i++)
//                {
//                        if (f_left_h[i] < f_min)
//                                continue;
//                        if (f_right_h[i] > f_max)
//                                continue;
//                        else
//                                out<<mu_hm_out[i]<<' ';
//                }
//
//        for(size_t i=0; i<mu_mp_out.size(); i++)
//                {
//                        if (f_left_m[i] < f_min)
//                                continue;
//                        if (f_right_m[i] > f_max)
//                                continue;
//                        else
//                                out<<mu_mp_out[i]<<' ';
//                }
//
//        for(size_t i=0; i<mu_mm_out.size(); i++)
//                {
//                        if (f_left_m[i] < f_min)
//                                continue;
//                        if (f_right_m[i] > f_max)
//                                continue;
//                        else
//                                out<<mu_mm_out[i]<<' ';
//                }
//

	for(size_t i=0; i<counts_hp.size(); i++)
        	{
	                if (f_left_h[i] < f_min)
                                continue;
                        if (f_right_h[i] > f_max)
                                continue;
			else
				out<<counts_hp[i]<<' ';
		}

//        for(size_t i=0; i<counts_hm.size(); i++)
//		{
//                        if (f_left_h[i] < f_min)
//                                continue;
//                        if (f_right_h[i] > f_max)
//                                continue;
//			else
//		                out<<counts_hm[i]<<' ';
//		}
//
//        for(size_t i=0; i<counts_mp.size(); i++)
//		{
//                        if (f_left_m[i] < f_min)
//                                continue;
//                        if (f_right_m[i] > f_max)
//                                continue;
//                        else
//		                out<<counts_mp[i]<<' ';
// 		}
//
//        for(size_t i=0; i<counts_mm.size(); i++)
//		{
//                        if (f_left_m[i] < f_min)
//                                continue;
//                        if (f_right_m[i] > f_max)
//                                continue;
//                        else
//		                out<<counts_mm[i]<<' ';  
//		}
}

string MyModel::description() const
{
	return string("");
}

