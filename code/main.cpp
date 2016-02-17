#include <iostream>
#include "DNest4/code/Start.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	CommandLineOptions options(argc, argv);
	Data::get_instance().load(options.get_data_file().c_str());
	Sampler<MyModel> sampler = setup<MyModel>(options);
	sampler.run();
	return 0;
}

