# NMRProbeCalib
Position calibration for NMR field probes. Written in C++ with a Matlab wrapper.
The mex file can be compiled using Matlab.
Uses a Newton Optimisation method to minimise the difference between the measured
probe frequencies and the frequencies estimated by the gradient fields.
Gradient fields up to 4th degree decomposition can be used.

See the ProgrammersManual.pdf for more information

Please cite:
Chang P, Eschelbach M, Syha R, Scheffler K and Henning A. "Impact of Gradient Nonlinearity on the Accuracy of NMR Field Camera Readouts"
Proc. of ISMRM 2015, Toronto, Canada.
