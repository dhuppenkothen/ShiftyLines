#include "Data.h"
#include "MyConditionalPrior.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstring>


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <CCfits>

using namespace std;


Data Data::instance;

Data::Data()
{

}

void Data::load_data(const char* datadir, const char* filename)
{
  pha_heg_p = load_fits(datadir, "1weak_4strong_heg_p1_src.pha");
  pha_meg_p = load_fits(datadir, "1weak_4strong_meg_p1_src.pha");
  pha_heg_m = load_fits(datadir, "1weak_4strong_heg_m1_src.pha");
  pha_meg_m = load_fits(datadir, "1weak_4strong_meg_m1_src.pha");

}

PHAData Data::load_fits(const char* datadir, const char* filename)
{
  char whole_file[256];
  strcpy(whole_file, datadir);
  strcat(whole_file, filename);
  strcpy(pha.filename, whole_file);

  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::unique_ptr<CCfits::FITS> input_file(new CCfits::FITS(pha.filename,CCfits::Read));

  // point to correct HDU extension
  CCfits::ExtHDU& spectrum = input_file->extension("SPECTRUM");

  // get out the data
  CCfits::Column& column = spectrum.column("CHANNEL");
  column.read(pha.channel, 1, column.rows());

  CCfits::Column& column2 = spectrum.column("COUNTS");
  column2.read(pha.counts, 1, column2.rows());

//  CCfits::Column& column3 = spectrum.column("BIN_LO");
//  column3.read(pha.bin_lo, 1, column3.rows());
 
//  CCfits::Column& column4 = spectrum.column("BIN_HI");
//  column4.read(pha.bin_hi, 1, column4.rows());
//
//  vector<double> bin_mid(pha.bin_lo.size(), 0.0);
//
//  // compute the middle of the energy bins
//  for(size_t i=0; i<pha.bin_lo.size(); i++){
//     bin_mid[i] = pha.bin_lo[i] + (pha.bin_hi[i] - pha.bin_lo[i])/2.0;
//  }
//
//  pha.bin_mid = bin_mid;
//
  string respfile, ancrfile;
  // need to read some keys, too!
  spectrum.readKey("RESPFILE", respfile);
  spectrum.readKey("ANCRFILE", ancrfile);

  string eunit_lo, eunit_hi;
 
  spectrum.readKey("TUNIT3", eunit_lo);
  spectrum.readKey("TUNIT4", eunit_hi);

  //cout<<"# Unit of the energy bins: "<<eunit_lo<<"."<<endl;

  pha.bin_lo_unit = "keV";//eunit_lo;
  pha.bin_hi_unit = "keV";//eunit_hi;

  pha.respfile = respfile;
  pha.ancrfile = ancrfile;

  cout<<"# Redistribution Matrix File: "<< pha.respfile <<"."<<endl;
  cout<<"# Ancillary Response File: "<< pha.ancrfile <<"."<<endl;

  // the number of counts
  pha.ncounts = pha.counts.size();

  // Print out the first value to check
  cout<<"# Found "<< pha.ncounts <<" data points in "<<pha.filename<<"."<<endl;

  // need to convert respfile from string to char
  const char *respfilechars = pha.respfile.c_str();

  // load associated RMF file
  pha.rmf = load_rmf(datadir, respfilechars); 

  // need to convert ancrfile from string to char
  const char *ancrfilechars = pha.ancrfile.c_str();

  pha.arf = load_arf(datadir, ancrfilechars);

  pha.bin_lo.assign(pha.counts.size(), 0.0);
  pha.bin_hi.assign(pha.counts.size(), 0.0);

  for (int i=0; i<pha.counts.size(); i++)
    {
    pha.bin_lo[i] = pha.arf.energ_lo[i];
    pha.bin_hi[i] = pha.arf.energ_hi[i];
    }

  return pha;

}



RMFData Data::load_rmf(const char* datadir, const char* filename)
{

  char whole_file[256];
  strcpy(whole_file, datadir);
  strcat(whole_file, filename); 
  strcpy(rmf.filename, whole_file);


  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::unique_ptr<CCfits::FITS> input_file(new CCfits::FITS(rmf.filename,CCfits::Read));

  // instantiate objects for temporary variables
  std::vector<std::valarray<int>> f_chan, n_chan;
  std::vector<std::valarray<double>>  mtx; 

  // point to correct HDU extension
  CCfits::ExtHDU& matrix = input_file->extension("MATRIX");
 
  // get out the data
  CCfits::Column& column = matrix.column("ENERG_LO");
  column.read(rmf.energ_lo, 1, column.rows());

  CCfits::Column& column2 = matrix.column("ENERG_HI");
  column2.read(rmf.energ_hi, 1, column2.rows());

  CCfits::Column& column3 = matrix.column("N_GRP");
  column3.read(rmf.n_grp, 1, column3.rows());

  CCfits::Column& column4 = matrix.column("F_CHAN");
  column4.readArrays(f_chan, 1, column4.rows());
 
  CCfits::Column& column5 = matrix.column("N_CHAN");
  column5.readArrays(n_chan, 1, column5.rows());

  // NOTE: Alternative names could be "SPECRESP MATRIX", 
  // "AXAF_RMF" and it could also be "HDUCLAS1" or "HDUCLAS2" 
  // with a key being "RESPONSE" or "RSP_MATRIX".
  // need to add that functionality some time.
  CCfits::Column& column_m = matrix.column("MATRIX");
  column_m.readArrays(mtx, 1, column_m.rows());

  cout<<"There are "<<rmf.energ_lo.size()<<" energy bins."<<endl;

  // somem keywords
  int detchans, tlmin;
  matrix.readKey("DETCHANS", detchans);
  try {
      // NOTE: Keyword hardcoded in, but probably shouldn't be!
      // In practice, need to figure out the number after TLMIN from 
      // FITS file. See Sherpa source code + notebook for details
      matrix.readKey("TLMIN4", tlmin);
      rmf.tlmin = tlmin;
      }
  catch(CCfits::HDU::NoSuchKeyword)
      { rmf.tlmin = 0; }

  rmf.detchans = detchans;

  rmf.offset = tlmin; 

  // point to correct HDU extension
  CCfits::ExtHDU& ebounds = input_file->extension("EBOUNDS");

  CCfits::Column& column6 = ebounds.column("E_MIN");
  column6.read(rmf.e_min, 1, column6.rows());

  CCfits::Column& column7 = ebounds.column("E_MAX");
  column7.read(rmf.e_max, 1, column7.rows());


  // do some magic on the arrays
  for(size_t i=0;i<rmf.n_grp.size();i++) 
    {
    if(rmf.n_grp[i]>0)
         {
         if (f_chan[i].size() > 1)
           {
           for(size_t j=0;j< f_chan[i].size(); j++)
              rmf.f_chan.push_back(f_chan[i][j]);
           
           }
        else
             rmf.f_chan.push_back(f_chan[i][0]);
       
        if (n_chan[i].size() > 1)
           {
           for(size_t j=0;j< n_chan[i].size(); j++)
              rmf.n_chan.push_back(n_chan[i][j]);
           }
        else
            rmf.n_chan.push_back(n_chan[i][0]);
       
        for(size_t j=0; j<mtx[i].size(); j++)
           rmf.matrix.push_back(mtx[i][j]);
  }
  }

  cout<<"There are "<<rmf.matrix.size()<<" elements in the response matrix."<<endl;

 
  return rmf; 
}

ARFData Data::load_arf(const char* datadir, const char* filename)
{

  char whole_file[256];
  strcpy(whole_file, datadir);
  strcat(whole_file, filename);
  strcpy(arf.filename, whole_file);

  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::unique_ptr<CCfits::FITS> input_file(new CCfits::FITS(arf.filename,CCfits::Read));

  // point to correct HDU extension
  // NOTE: This could also be "AXAF_ARF"
  CCfits::ExtHDU& specresp = input_file->extension("SPECRESP");

  // get out the data
  CCfits::Column& column = specresp.column("ENERG_LO");
  column.read(arf.energ_lo, 1, column.rows());

  CCfits::Column& column2 = specresp.column("ENERG_HI");
  column2.read(arf.energ_hi, 1, column2.rows());

  CCfits::Column& column3 = specresp.column("SPECRESP");
  column3.read(arf.specresp, 1, column3.rows());

  try
    {
    CCfits::Column& column4 = specresp.column("BIN_LO");
    column4.read(arf.bin_lo, 1, column4.rows());
  
    CCfits::Column& column5 = specresp.column("BIN_HI");
    column5.read(arf.bin_hi, 1, column5.rows());
    }
  catch(CCfits::HDU::NoSuchKeyword)
    { 
    cout<<"No keywords BIN_LO and BIN_HI"<<endl;
    }

  cout<<"There are "<<arf.specresp.size()<<" elements in the ARF."<<endl;
  
  return arf;
} 

void Data::load_lines(const char* filename)
{
	double conv = 12.3984191;

	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		exit(1);
	}

	lines.clear();

	// ASSUME DATA IS IN KEV. IF NOT, THEN CONVERT IT TO KEV!
        double temp1;
        while(fin>>temp1)
                //if (pha.bin_lo_unit == "angstrom")
	        //        lines.push_back(conv/temp1);
		//else if (pha.bin_lo_unit == "kev")
			lines.push_back(temp1);
		//else
		//	throw CCfits::Column::InvalidDataType(); 

	nlines = lines.size();

        fin.close();
        cout<<"# Found "<<lines.size()<<" line positions in "<<filename<<"."<<endl;
}


void Data::load(const char* filename)
{

        double conv = 12.3984191;

	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		exit(1);
	}

	f_left.clear();
	f_right.clear();
	y.clear();
	yerr.clear();

	double temp1, temp2, temp3, temp4;
	while(fin>>temp1 && fin>>temp2 && fin>>temp3 && fin>>temp4)
	{
		f_left.push_back(temp1);
		f_right.push_back(temp2);
		y.push_back(temp3);
		yerr.push_back(temp4);
	}

	fin.close();
	cout<<"# Found "<<f_left.size()<<" points in file "<<filename<<"."<<endl;

	compute_summaries();
}


void Data::compute_summaries()
{
	f_min = 1.83;
	f_max = 2.03;
	//f_min = 1.746256187217072;
	//f_max = 2.0493255531922747;

	// ASSUME DATA IS IN KEV!!
	
	//const double l_min = lines[0];
	//const double l_max = lines[lines.size()-1];

	//const double dmin = -0.01;
	//const double dmax = 0.01;

	//f_min = l_min/(1. + dmax) - 0.01;
	//f_max = l_max/(1 + dmin) - 0.01;
	f_range = f_max - f_min;

	cout<<"Total energy range covered: "<<f_range<<"."<<endl;

	// Left and right edges of the data bins
	f_mid.assign(f_left.size(), 0.);
	df.assign(f_left.size(), 0.);
	for(size_t i=0; i<f_left.size(); i++)
	{
		df[i] = f_right[i] - f_left[i];
		f_mid[i] = f_left[i] + 0.5*df[i];
	}

	min_df = *min_element(df.begin(), df.end());
}


