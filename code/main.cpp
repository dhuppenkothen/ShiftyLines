#include <iostream>
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Data.h"
#include <string>
using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	// parse command line options
	CommandLineOptions options(argc, argv);
 
        //char datadir[50];
        //strcpy(datadir, "../data/cyg_daniela/");
	std::string datadir = "../data/cygx1_sims/";

	// get the data filename from the command line
	Data::get_instance().load_data(datadir.c_str(), options.get_data_file().c_str());

        cout<<"Loaded FITS files with data ..."<<endl;

	// I need the old text file to set min/max energy boundaries for the spectrum
        char text_file[128];
        strcpy(text_file, "../data/test_noshift1.txt");

        Data::get_instance().load(text_file);
//	 Data::get_instance().load(options.get_data_file().c_str());

        cout<<"Loaded text file data ..."<<endl;


	// file with line positions in same units as data; CURRENTLY HARDCODED!
	Data::get_instance().load_lines("../data/si_lines_kev.txt");
 
//        cout<<"Loaded line data ..."<<endl;
 
	// sample!
	Sampler<MyModel> sampler = setup<MyModel>(options);
  
//        cout<<"Set up sampler ..."<<endl;

	sampler.run();
	return 0;
}

