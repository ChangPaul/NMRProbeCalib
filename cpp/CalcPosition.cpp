/*
 * CalcPosition.cpp - Calculates the probe positions from the measured
 *					  frequencies (normalised).
 *					  Uses a Newton optimisation method
 * 
 * MEX-file for MATLAB.
 *
 */

#include "mex.h"
#include "mex_err.h"
#include "OptimalPosition/OptimalPositions.h"
#include <vector>

using namespace std;

//==================================================================
//                       Function Implementations
//==================================================================

/*
 * FUNCTION:    CalcPosition
 * DESCRIPTION: Calculates the optimal positions of the probes
 *              by calling the OptimalPositions class.
 *              The function mainly organises the variables so that
 *              the OptimalPositions class can use them for the
 *              calculations.
 * INPUTS:      opt_pos    - linear array of the optimal x-, y-,
 *                           z-positions of the probes
 *              initpos    - linear array of the initial x-, y-,
 *                           z-positions of the probes
 *              nprobe     - number of probes
 *              gradcoeffs - linear array of the grad coeffs
 *              numsh      - number of sph harm coeffs
 */
void CalcPosition (double* opt_pos, double* initpos, int nprobe,
                   double* gradcoeffs, int numsh)
{
	// Read gradient coefficients
	fvec xgrad, ygrad, zgrad;
	if (numsh == 0 || gradcoeffs == NULL)
	{
		// Default gradient coefficients
		xgrad.resize(3);	xgrad[0] = 1.0;
		ygrad.resize(3);	ygrad[1] = 1.0;
		zgrad.resize(3);	zgrad[2] = 1.0;
	}
	else
	{
		// Store grad coefficients from input array
		for (int i = 0; i != numsh; ++i)
		{
			xgrad.push_back (gradcoeffs[i*3 + 0]);
			ygrad.push_back (gradcoeffs[i*3 + 1]);
			zgrad.push_back (gradcoeffs[i*3 + 2]);
		}
	}
	
	// Instantiate an "OptimalPositions" object with the given
	// gradient coefficients.
	fvec3d freq;
	OptimalPositions optpos (xgrad, ygrad, zgrad);
	
	// Calculate the optimal positions for each probe.
	for (int i = 0; i != nprobe; ++i)
	{
		freq[0] = initpos [i*3 + 0];
		freq[1] = initpos [i*3 + 1];
		freq[2] = initpos [i*3 + 2];
		optpos.Calculate (freq);
		opt_pos [i*3 + 0] = freq [0];
		opt_pos [i*3 + 1] = freq [1];
		opt_pos [i*3 + 2] = freq [2];
	}
}

//==================================================================
//                           Main Function
//==================================================================

/*
 * Calculates the positions of the probes using a Newton's optimisation
 * algorithm to find the best fit for the frequencies.
 * This function mainly checks the input arguments and reorganises them
 * for the CalcPosition function (above).
 */
 
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Declare and initialise variables
	double *pfreq, *gradcoeffs = NULL;
	double *pos;
	mwSize nprobe, numsh;
	int err = 0;

	// Obtain variables from the input arguments
	// Freq data is stored in the first element of the arg array.
	// Grad coeffs are stored in the second element of the arg array.
	switch (nrhs)
	{
	case 2: // Gradient coefficients
			if (!mxIsDouble (prhs[1]) || mxIsComplex(prhs[1])) 	err = 0x102;
			else if (mxGetM (prhs[1]) != 3) 					err = 0x103;
			else gradcoeffs = mxGetPr (prhs[1]);
	case 1: // Probe frequency data (Nx3 matrix)
			if (!mxIsDouble (prhs[0]) || mxIsComplex(prhs[0])) 	err = 0x102;
			else if (mxGetM (prhs[0]) != 3) 					err = 0x103;
			else pfreq = mxGetPr (prhs[0]);
			break;
	case 0: err = 0x100; break;
	default: err = 0x101;
	}
	
	// Check if output variable was assigned
	if (nlhs != 1) err = 0x200;
	error (err, "CalcPosition");
	
	// Extract the number of probes and the number of sph harm coeffs
	nprobe = (mwSize) mxGetN (prhs[0]);
	numsh  = (nrhs > 1) ? (mwSize) mxGetN (prhs[1]) : 0;
	
	// Create output
	plhs[0] = mxCreateDoubleMatrix (3, nprobe, mxREAL);
	pos = mxGetPr (plhs[0]);
	
	// Read probe frequencies and calculate optimum positions
	CalcPosition (pos, pfreq, nprobe, gradcoeffs, numsh);
}
