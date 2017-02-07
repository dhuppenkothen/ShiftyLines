import shutil
import argparse
def move_dnest_output(froot, dnest_dir="./"):
    shutil.move(dnest_dir+"posterior_sample.txt", froot+"_posterior_sample.txt")
    shutil.move(dnest_dir+"sample.txt", froot+"_sample.txt")
    shutil.move(dnest_dir+"sample_info.txt", froot+"_sample_info.txt")
    shutil.move(dnest_dir+"weights.txt", froot+"_weights.txt")
    shutil.move(dnest_dir+"levels.txt", froot+"_levels.txt")
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--froot", action="store", required=True, dest="froot", type=str)
    parser.add_argument("-d", "--dnest_dir", action="store", dest="dnest_dir", required=False, type=str, default="./")
    clargs = parser.parse_args()

    move_dnest_output(clargs.froot, clargs.dnest_dir) 
