#include <iostream>
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	// parse command line options
	CommandLineOptions options(argc, argv);

	// get the data filename from the command line
	Data::get_instance().load(options.get_data_file().c_str());

	// file with line positions in same units as data; CURRENTLY HARDCODED!
	Data::get_instance().load_lines("../data/si_lines.txt");
 
	// sample!
	Sampler<MyModel> sampler = setup<MyModel>(options);
	sampler.run();
	return 0;
}

