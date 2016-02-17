import numpy as np
import matplotlib.pyplot as plt

import os
import argparse

#data = loadtxt('../data/test_data3.dat')
#posterior_sample = atleast_2d(loadtxt('posterior_sample.txt'))

def display_samples(filename, posterior_dir, save_frames=False):
    """
    Utility function for displaying samples from the posterior along 
    with the data. 

    Parameters
    ----------
    filename: str
        full path of the file containing the data

    posterior_dir: str
        The directory containing the posterior_samples.txt file

    save_frames: bool
        Flag determining whether to save the frames being plotted 
        to files. Default: False


    """
    
    data = np.loadtxt(filename)
    posterior_sample = np.atleast_2d(np.loadtxt(posterior_dir+"posterior_sample.txt"))

    if save_frames:
        os.system('rm '+ posterior_dir + 'Frames/*.png')
        try:
            os.mkdir(posterior_dir+"Frames/")        
        except OSError:	
            print("Directory already exists!")

    plt.ion()
    for i in range(0, posterior_sample.shape[0]):
        plt.hold(False)
        plt.plot(data[:,0], data[:,2], linestyle="steps-mid")
        plt.hold(True)
        plt.plot(data[:,0], posterior_sample[i, -data.shape[0]:], 'r')
        plt.xlabel('Wavelength [Angstrom]', fontsize=16)
        plt.ylabel('Flux', fontsize=16)
        plt.draw()
        if save_frames:
            plt.savefig(posterior_dir+'Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
            print(posterior_dir+'Frames/' + '%0.4d'%(i+1) + '.png')

    plt.ioff()
    plt.show()
    return
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="Display data and posterior samples.")


    parser.add_argument('-f', '--filename', action="store", dest="filename",
                        required=True, help="Filename of the data file to be plotted.")
    parser.add_argument('-p', '--posterior_dir', action="store", dest="posterior_dir", required=False, default="./",
                        help="Directory with the posterior_samples.txt file.")

    parser.add_argument('-s', '--save', action="store", dest="save_frames", type=bool, required=False, default=False, 
                        help="Optional boolean flag determining whether to save each sample plot to disk.")


    clargs = parser.parse_args()

    filename = clargs.filename
    posterior_dir = clargs.posterior_dir
    save_frames = clargs.save_frames 

    display_samples(filename, posterior_dir, save_frames)
