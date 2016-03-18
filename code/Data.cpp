#include "Data.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <CCfits>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load_data(const char* filename)
{


  CCfits::FITS::setVerboseMode(true);

  // Open file for reading
  std::auto_ptr<CCfits::FITS> input_file(new CCfits::FITS(filename,CCfits::Read));

  // point to correct HDU extension
  CCfits::ExtHDU& spectrum = input_file->extension("SPECTRUM");

  // get out the data
  CCfits::Column& column = spectrum.column("CHANNEL");
  column.read(channel, 1, column.rows());
  
  CCfits::Column& column2 = spectrum.column("COUNTS");
  column.read(counts, 1, column2.rows());

  CCfits::Column& column3 = spectrum.column("BIN_LO");
  column.read(bin_lo, 1, column3.rows());

  CCfits::Column& column4 = spectrum.column("BIN_HI");
  column.read(bin_hi, 1, column4.rows());

  // need to read some keys, too!

  
  // the number of counts
  ncounts = counts.size();

  // Print out the first value to check
  cout<<"# Found "<< ncounts <<" line positions in "<<filename<<"."<<endl;


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


