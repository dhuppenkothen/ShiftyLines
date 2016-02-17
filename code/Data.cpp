#include "Data.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
using namespace std;

Data Data::instance;

Data::Data()
:line_pos {.64794778337, 6.18222706825,
                    6.71993761152, 6.78432363818,
                    6.85605388882, 6.9281555301,
                    6.99681113087, 7.05787012547}
,nlines(line_pos.size())
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Failed to open file "<<filename<<"."<<endl;
		return;
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

}

