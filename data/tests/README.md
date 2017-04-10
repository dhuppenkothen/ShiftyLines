Testing the Model
=================

This folder contains some results of models I used to test whether the algorithm does what it's supposed to.
The files are named after the branches of ShiftyLines that produced this output.

* test1 : only `HEG_P` data, model includes OU process, prints out all model components
* test2: same as above, but the lines are added to the log-background, such that they always must be positive
* test3: without OU process, so I can check some intuition I have about the RMF smoothing out the model and the OU process picking up some slope it maybe shouldn't.
