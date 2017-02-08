#ifndef _Data_
#define _Data_

#include <vector>
#include <string> 
#include <cstring>
#include <valarray>


typedef struct {
   // keep file name for later use
   char  filename[200];

   // the vectors with some of the data stuff;
   std::vector<double> energ_lo, energ_hi;
 
   std::vector<int> f_chan, n_chan, n_grp;
   std::vector<double>  matrix;
   //std::vector<std::valarray<float>> f_chan, n_chan, matrix;

   std::vector<int> e_min, e_max;

   int detchans, tlmin, offset;

}RMFData;

typedef struct {

  // store file name just in case
  char filename[200];
  
  // keyword
  std::string exposure;

  // the spectral response, energy bins
  std::vector<double> energ_lo, energ_hi, specresp, bin_lo, bin_hi;
  

}ARFData;
 
typedef struct {
   // keep file name for later use
   char  filename[200];

   // the data: channels, counts, energy bin edges
   std::vector<double> channel, counts, bin_lo, bin_mid, bin_hi;
   
   // the exposure and the number of counts
   double exposure, ncounts;

   // file names of ARF and RMF files
   std::string ancrfile, respfile;

   // units of bin_lo and bin_hi:
   std::string bin_lo_unit, bin_hi_unit;

   // Data structure with the RMF
   RMFData rmf;
   
   // Data structure with the ARF;
   ARFData arf;

}PHAData;


class Data
{
	private:
		// defining the data coming out of the fits files
		//double exposure, ncounts;
		//std::string ancrfile, respfile;

		//std::vector<double> channel, counts, bin_lo, bin_hi;
	
		PHAData pha;
   		RMFData rmf;
		ARFData arf;	

		std::vector<double> f_left, f_right, y, yerr;
		std::vector<double> f_mid, df;

		std::vector<double> lines;
		int nlines;

		// Some useful summaries
		double f_min, f_max, f_range, min_df;
		void compute_summaries();
	public:
		Data();
		void load(const char* filename);
		void load_lines(const char* filename);
      		void load_data(const char* datadir, const char* filename);
		RMFData load_rmf(const char* datadir, const char* filename);
                ARFData load_arf(const char* datadir, const char* filename);
                PHAData load_fits(const char* datadir, const char* filename);
		// Getters
		const std::vector<double>& get_f_mid() const { return f_mid; }
		const std::vector<double>& get_df() const { return df; }
		const std::vector<double>& get_f_left() const { return f_left; }
		const std::vector<double>& get_f_right() const { return f_right; }
		const std::vector<double>& get_y() const { return y; }
		const std::vector<double>& get_yerr() const { return yerr; }
		const std::vector<double>& get_line_pos() const { return lines; }
                const PHAData& get_pha() const { return pha; }

		double get_f_min() const { return f_min; }
		double get_f_max() const { return f_max; }
		double get_f_range() const { return f_range; }
		double get_min_df() const { return min_df; }
		const int& get_nlines() const { return nlines; }
		
	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

