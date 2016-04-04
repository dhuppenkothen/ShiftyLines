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
 
//        char datadir[50];
//        strcpy(datadir, "../data/cyg_daniela/");

	// get the data filename from the command line
//	Data::get_instance().load_data(datadir, options.get_data_file().c_str());
//        cout<<"Loaded FITS file data ..."<<endl;

//        char text_file[128];
//        strcpy(text_file, "../data/test_noshift1.txt");

//        Data::get_instance().load(text_file);
	 Data::get_instance().load(options.get_data_file().c_str());

        cout<<"Loaded text file data ..."<<endl;


	// file with line positions in same units as data; CURRENTLY HARDCODED!
	Data::get_instance().load_lines("../data/si_lines.txt");
 
//        cout<<"Loaded line data ..."<<endl;
 
	// sample!
	Sampler<MyModel> sampler = setup<MyModel>(options);
  
//        cout<<"Set up sampler ..."<<endl;

	sampler.run();
	return 0;
}

