import dnest4
import display
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description="Display diagnostics and results of the sampling.")

parser.add_argument('-f', '--filename', action="store", dest="filename",
                        required=True, help="Filename of the data file to be plotted.")

parser.add_argument('-p', '--posterior_dir', action="store", dest="posterior_dir", required=False, default="./",
                        help="Directory with the posterior_samples.txt file.")

clargs = parser.parse_args()

filename = clargs.filename
posterior_dir = clargs.posterior_dir

dnest4.postprocess()
display.display_samples(filename, posterior_dir)

