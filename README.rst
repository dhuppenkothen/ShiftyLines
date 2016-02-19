Shifty Lines
============

A Bayesian hierarchical inference model for X-ray spectra with 
emission and/or absorption lines that may be shifted by some Doppler 
shift. Infers the number of Doppler shifts present in the data along 
with line intensities and widths. 

Uses DNest4 and RJObject for reversible-jump MCMC.

Prerequisites
-------------
* gcc
* python3
* numpy
* `DNest4 <https://github.com/eggplantbren/DNest4>`_

Installing and Running the Model
--------------------

First, make sure that `DNest4 <https://github.com/eggplantbren/DNest4>`_ is installed.
Also, make sure you have a variable `DNEST4_PATH` set to the directory that 
the DNest4 repository lives in, and have the `DNest4/code/` directory in your 
`PYTHONPATH` variable.

You should now be able to go into `ShiftyLines/code/` and type `make` to compile 
the model.

Run  the model with 

    >>./main -d path/to/data/file -t 1

Use `-d` to set the path to the name of the file to be modelled. `-t` sets the 
number of parallel threads to be used: use `-t 1` for a single thread, or 
more for parallel threads. For best results, choose something like `-t 8` or 
something similar.


Copyright
---------

All content © 2016 the authors. The code is distributed under the MIT license.

Pull requests are welcome! If you are interested in the further development of
this project, please `get in touch via the issues
<https://github.com/dhuppenkothen/ShiftyLines/issues>`_!

