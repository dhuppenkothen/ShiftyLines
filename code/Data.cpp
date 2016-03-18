#include "Data.h"
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
  char whole_file[256];
  strcpy(whole_file, datadir);
  strcat(whole_file, filename);
  strcpy(pha.filename, whole_file);

  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::auto_ptr<CCfits::FITS> input_file(new CCfits::FITS(pha.filename,CCfits::Read));

  // point to correct HDU extension
  CCfits::ExtHDU& spectrum = input_file->extension("SPECTRUM");

  // get out the data
  CCfits::Column& column = spectrum.column("CHANNEL");
  column.read(pha.channel, 1, column.rows());
  
  CCfits::Column& column2 = spectrum.column("COUNTS");
  column.read(pha.counts, 1, column2.rows());

  CCfits::Column& column3 = spectrum.column("BIN_LO");
  column.read(pha.bin_lo, 1, column3.rows());

  CCfits::Column& column4 = spectrum.column("BIN_HI");
  column.read(pha.bin_hi, 1, column4.rows());

  vector<double> bin_mid(pha.bin_lo.size(), 0.0);

  // compute the middle of the energy bins
  for(size_t i=0; i<pha.bin_lo.size(); i++){
     bin_mid[i] = pha.bin_hi[i] - pha.bin_lo[i];   
  }

  pha.bin_mid = bin_mid;

  string respfile, ancrfile;
  // need to read some keys, too!
  spectrum.readKey("RESPFILE", respfile);
  spectrum.readKey("ANCRFILE", ancrfile);

  string eunit_lo, eunit_hi;
  spectrum.readKey("TUNIT2", eunit_lo);
  spectrum.readKey("TUNIT3", eunit_hi);

  cout<<"# Unit of the energy bins: "<<eunit_lo<<"."<<endl;

  pha.bin_lo_unit = eunit_lo;
  pha.bin_hi_unit = eunit_hi;

  pha.respfile = respfile;
  pha.ancrfile = ancrfile;

  cout<<"# Redistribution Matrix File: "<< pha.respfile <<"."<<endl;
  cout<<"# Ancillary Response File: "<< pha.ancrfile <<"."<<endl;

  // the number of counts
  pha.ncounts = pha.counts.size();

  // Print out the first value to check
  cout<<"# Found "<< pha.ncounts <<" line positions in "<<pha.filename<<"."<<endl;

  const char *respfilechars = pha.respfile.c_str();

  rmf = load_rmf(datadir, respfilechars); 
  // load the associated rmf files

}

RMFData Data::load_rmf(const char* datadir, const char* filename)
{

  char whole_file[256];
  strcpy(whole_file, datadir);
  strcat(whole_file, filename); 
  strcpy(rmf.filename, whole_file);

  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::auto_ptr<CCfits::FITS> input_file(new CCfits::FITS(rmf.filename,CCfits::Read));

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
  column4.read(rmf.f_chan, 1, column4.rows());


  //CCfits::Column& column5 = matrix.column("N_CHAN");
  //column5.read(rmf.n_chan, 1, column5.rows());

  cout<<"There are "<<rmf.energ_lo.size()<<" energy bins."<<endl;

  // somem keywords
  int detchans, tlmin;
  matrix.readKey("DETCHANS", detchans);
  matrix.readKey("TLMIN", tlmin);
  
  rmf.detchans = detchans;
  rmf.tlmin = tlmin;

  rmf.offset = tlmin; 

  // point to correct HDU extension
  CCfits::ExtHDU& ebounds = input_file->extension("MATRIX");

  CCfits::Column& column6 = ebounds.column("E_MIN");
  column6.read(rmf.e_min, 1, column6.rows());

  CCfits::Column& column7 = ebounds.column("E_MAX");
  column7.read(rmf.e_max, 1, column7.rows());
  
  return rmf; 
}

void Data::load_lines(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		exit(1);
	}

	lines.clear();

        double temp1;
        while(fin>>temp1)
                lines.push_back(temp1);

	nlines = lines.size();

        fin.close();
        cout<<"# Found "<<lines.size()<<" line positions in "<<filename<<"."<<endl;
	
	

}



void Data::load(const char* filename)
{
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
	f_min = *min_element(f_left.begin(), f_left.end());
	f_max = *max_element(f_right.begin(), f_right.end());
	f_range = f_max - f_min;

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


