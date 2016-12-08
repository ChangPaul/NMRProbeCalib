/*
 * OptimalPositions.cpp
 *
 *  Created on: Jun 11, 2015
 *      Author: changp
 */

#include "mex.h"
#include "OptimalPositions.h"

//===============================================================
// 						PUBLIC METHODS
//===============================================================

/*
 * FUNCTION:	(constructor)
 * INPUTS:		xgrad - sph harmonic coeffs of x-gradient
 * 				ygrad - sph harmonic coeffs of y-gradient
 * 				zgrad - sph harmonic coeffs of z-gradient
 */
OptimalPositions::OptimalPositions( const fvec &xgrad, const fvec &ygrad, const fvec &zgrad )
{
	num_fn = 25;

	kx = xgrad;	kx.resize( num_fn );
	ky = ygrad;	ky.resize( num_fn );
	kz = zgrad;	kz.resize( num_fn );
}

/*
 * FUNCTION:	Calculate
 * DESCRIPTION:	Calculate the optimal position that matches the target
 * 				frequencies.
 * INPUTS:		init - initial position.
 * OUTPUTS:		init - optimal position (note that the variable is
 * 					   overwritten).
 * RETURNS:		Final error between estimate and target.
 */
float OptimalPositions::Calculate( fvec3d init )
{
	fvec3d freq = { init[0], init[1], init[2] };
	fvec3d pos  = { init[0], init[1], init[2] };
	fvec3d step;

	// Iterate until convergence or timeout (max 100 iterations)
	int itr = 0, timeout = 20;
	float prev = 1000.0, err = 100.0;
	float threshold = 1e-6;
	while( abs( prev - err ) > threshold && ++itr < timeout )
	{
		prev = err;
		err = OptimStep( pos, freq, step );
		pos[0] += step[0];
		pos[1] += step[1];
		pos[2] += step[2];
	}

	// Return output
	init[0] = pos[0];
	init[1] = pos[1];
	init[2] = pos[2];
	
	// Return final error
	return Fval( pos, freq );
}

//===============================================================
// 						PRIVATE METHODS
//===============================================================

/*
 * FUNCTION:	OptimStep
 * DESCRIPTION:	Calculates a step direction for the optimisation
 * 				algorithm.
 * INPUTS:		pos  - current position.
 * 				freq - measured freq (i.e. target values).
 * OUTPUTS:		step - step direction.
 * RETURNS:		Either conditioning of the Hessian (for Newton method)
 * 				or the scaling of the step (for gradient method using
 * 				Armijo line search method).
 */
float OptimalPositions::OptimStep( const fvec3d &pos, const fvec3d &freq, fvec3d step )
{
	float alpha;
	fvec3d dir;
	fmat3d hessian, invHess;

	dFval( pos, freq, dir );
	
// Newton's method
#if ALGORITHM == 0
	// Calculate inverse Hessian matrix
	ddFval( pos, freq, hessian );
	alpha = invertMatrix( hessian, invHess );

	// Calculate step
	step[0] = -dir[0]*invHess[0][0] - dir[1]*invHess[0][1] - dir[2]*invHess[0][2];
	step[1] = -dir[0]*invHess[1][0] - dir[1]*invHess[1][1] - dir[2]*invHess[1][2];
	step[2] = -dir[0]*invHess[2][0] - dir[1]*invHess[2][1] - dir[2]*invHess[2][2];
	
// Gradient descent method
#elif ALGORITHM == 1

	// Armijo line search
	step[0] = -dir[0]; step[1] = -dir[1]; step[2] = -dir[2];
	float a = 2.0, sigma = 0.8, gamma = 0.25;
	float err;
	fvec3d pos_next;
	
	float alpha = a/sigma;
	do
	{
		alpha *= sigma;
		pos_next[0] = pos[0] + alpha*step[0];
		pos_next[1] = pos[1] + alpha*step[1];
		pos_next[2] = pos[2] + alpha*step[2];
		err  = Fval( pos_next, freq ) - Fval( pos, freq );
		err -= gamma*alpha*( step[0]*dir[0] + step[1]*dir[1] + step[2]*dir[2] );
	} while( err > 0 );

	step[0] *= alpha; step[1] *= alpha; step[2] *= alpha;
#endif

	return alpha;
}

/*
 * FUNCTION:	EstimateError
 * DESCRIPTION:	Calculates the error between the estimated and
 * 				measured frequencies.
 * INPUTS:		pos  - current position.
 * 				freq - measured frequencies (i.e. target values).
 * OUTPUTS:		err  - difference between estimated and measured freq.
 */
void OptimalPositions::EstimateError( const fvec3d &pos, const fvec3d &freq, fvec3d err )
{
	float resultx = 0.0, resulty = 0.0, resultz = 0.0;
	for( int i = 0; i != num_fn; ++i )
	{
		float val = sphFn[i]( pos[0], pos[1], pos[2] );
		resultx += kx[i]*val;
		resulty += ky[i]*val;
		resultz += kz[i]*val;
	}
	err[0] = resultx - freq[0];
	err[1] = resulty - freq[1];
	err[2] = resultz - freq[2];
}

/*
 * FUNCTION:	Fval
 * DESCRIPTION:	Evaluates the objective function for the current
 * 				position:
 * 					|| est_freq - meas_freq ||^2
 * INPUTS:		pos  - position to be evaluated.
 * 				freq - measured frequency (i.e. target values).
 * RETURNS:		Evaluated objective function.
 */
float OptimalPositions::Fval( const fvec3d &pos, const fvec3d &freq )
{
	fvec3d err;
	EstimateError( pos, freq, err );
	return err[0]*err[0] + err[1]*err[1] + err[2]*err[2];
}

/*
 * FUNCTION:	dFval
 * DESCRIPTION:	Calculates the Jacobian of the objective function:
 * 					dF = [ dF/dx, dF/dy, dF/dz ]';
 * 				The objective function is the L2-norm squared of the
 * 				estimate error:
 * 					|| est_freq - meas_freq ||^2
 * INPUTS:		pos  - current position.
 * 				freq - measured frequencies (i.e. target values).
 * OUTPUTs:		dx   - the Jacobian vector.
 */
void OptimalPositions::dFval( const fvec3d &pos, const fvec3d &freq, fvec3d dF )
{
	fvec3d err;
	fvec3d dFx, dFy, dFz;

	EstimateError( pos, freq, err );
	dFdx( pos, kx, dFx );
	dFdx( pos, ky, dFy );
	dFdx( pos, kz, dFz );

	dF[0] = err[0]*dFx[0] + err[1]*dFy[0] + err[2]*dFz[0];
	dF[1] = err[0]*dFx[1] + err[1]*dFy[1] + err[2]*dFz[1];
	dF[2] = err[0]*dFx[2] + err[1]*dFy[2] + err[2]*dFz[2];
}

/*
 * FUNCTION:	ddFval
 * DESCRIPTION: Calculates the Hessian matrix of the given function at
 * 				the given position.
 * INPUTS:		pos  - position to be evaluated.
 * 				freq - freq measured at the target position.
 * OUTPUTS:		ddF  - Hessian matrix.
 */
void OptimalPositions::ddFval( const fvec3d &pos, const fvec3d &freq, fmat3d ddF )
{
	// Calculate Jacobians of the x-, y-, z-grad
	fvec3d dFx, dFy, dFz;
	dFdx( pos, kx, dFx );
	dFdx( pos, ky, dFy );
	dFdx( pos, kz, dFz );

	// Calculate the error between current estimate and actual freq
	fvec3d err;
	EstimateError( pos, freq, err );

	//-----------------------------------------------
	// Calculate the first term of the Hessian matrix
	ddF[0][0] = dFx[0]*dFx[0] + dFy[0]*dFy[0] + dFz[0]*dFz[0];
	ddF[0][1] = dFx[0]*dFx[1] + dFy[0]*dFy[1] + dFz[0]*dFz[1];
	ddF[0][2] = dFx[0]*dFx[2] + dFy[0]*dFy[2] + dFz[0]*dFz[2];
	ddF[1][0] = dFx[1]*dFx[0] + dFy[1]*dFy[0] + dFz[1]*dFz[0];
	ddF[1][1] = dFx[1]*dFx[1] + dFy[1]*dFy[1] + dFz[1]*dFz[1];
	ddF[1][2] = dFx[1]*dFx[2] + dFy[1]*dFy[2] + dFz[1]*dFz[2];
	ddF[2][0] = dFx[2]*dFx[0] + dFy[2]*dFy[0] + dFz[2]*dFz[0];
	ddF[2][1] = dFx[2]*dFx[1] + dFy[2]*dFy[1] + dFz[2]*dFz[1];
	ddF[2][2] = dFx[2]*dFx[2] + dFy[2]*dFy[2] + dFz[2]*dFz[2];

	//--------------------------------------------------
	// Calculate the second derivative for each gradient
	vector< float > kdx, kdy, kdz;
	vector< float > kdxdx, kdydx, kdzdx;
	vector< float > kdxdy, kdydy, kdzdy;
	vector< float > kdxdz, kdydz, kdzdz;

	// Iterate over x-, y-, z-gradients
	fmat3d secDeriv;
	for( int g = 0; g != 3; ++g )
	{
		switch( g )
		{
		case 0: derivCoeffs( kx, &kdx, &kdy, &kdz ); break;
		case 1: derivCoeffs( ky, &kdx, &kdy, &kdz ); break;
		case 2: derivCoeffs( kz, &kdx, &kdy, &kdz ); break;
		}

		// Calculate coefficients of the second derivatives
		derivCoeffs( kdx, &kdxdx, &kdxdy, &kdxdz );
		derivCoeffs( kdy, &kdydx, &kdydy, &kdydz );
		derivCoeffs( kdz, &kdzdx, &kdzdy, &kdzdz );

		// Evaluate second derivatives at the given position
		secDeriv[0][0] = 0.0;		secDeriv[0][1] = 0.0;		secDeriv[0][2] = 0.0;
		secDeriv[1][0] = 0.0;		secDeriv[1][1] = 0.0;		secDeriv[1][2] = 0.0;
		secDeriv[2][0] = 0.0;		secDeriv[2][1] = 0.0;		secDeriv[2][2] = 0.0;
		for( unsigned i = 0; i != kdx.size(); ++i )
		{
			float val = sphFn[i]( pos[0], pos[1], pos[2] );
			secDeriv[0][0] += kdxdx[i]*val;		secDeriv[1][0] += kdydx[i]*val;		secDeriv[2][0] += kdzdx[i]*val;
			secDeriv[0][1] += kdxdy[i]*val;		secDeriv[1][1] += kdydy[i]*val;		secDeriv[2][1] += kdzdy[i]*val;
			secDeriv[0][2] += kdxdz[i]*val;		secDeriv[1][2] += kdydz[i]*val;		secDeriv[2][2] += kdzdz[i]*val;
		}

		// Update the second term of the Hessian
		ddF[0][0] += secDeriv[0][0]*err[g];	ddF[0][1] += secDeriv[0][1]*err[g];	ddF[0][2] += secDeriv[0][2]*err[g];
		ddF[1][0] += secDeriv[1][0]*err[g];	ddF[1][1] += secDeriv[1][1]*err[g];	ddF[1][2] += secDeriv[1][2]*err[g];
		ddF[2][0] += secDeriv[2][0]*err[g];	ddF[2][1] += secDeriv[2][1]*err[g];	ddF[2][2] += secDeriv[2][2]*err[g];
	}
}

/*
 * FUNCTION:	dFdx
 * DESCRIPTION:	Calculates the value of the Jacobian at the given point.
 * INPUTS:		pos - position to be evaluated.
 * 				k   - sph harmonic coeffs of the function.
 * OUTPUTS:		dx  - Jacobian at the given position
 */
void OptimalPositions::dFdx( const fvec3d &pos, const vector< float > &k, fvec3d dx )
{
	// Calculate the coeffs of the deriv functions
	vector< float > kdx, kdy, kdz;
	derivCoeffs( k, &kdx, &kdy, &kdz );

	// Evaluate deriv functions at the given position
	float resultx = 0.0, resulty = 0.0, resultz = 0.0;
	for( unsigned i = 0; i != k.size(); ++i )
	{
		float val = sphFn[i]( pos[0], pos[1], pos[2] );
		resultx += kdx[i]*val;
		resulty += kdy[i]*val;
		resultz += kdz[i]*val;
	}

	dx[0] = resultx; dx[1] = resulty; dx[2] = resultz;
}

/*
 * FUNCTION:	derivCoeffs
 * DESCRIPTION:	Calculates the coefficients of the derivatives of the
 * 				spherical harmonic functions.
 * INPUTS:		k  - coeffs of sph harm function
 * OUTPUT:		kx - coeffs of x derivative.
 * 				ky - coeffs of y derivative.
 * 				kz - coeffs of z derivative.
 */
void OptimalPositions::derivCoeffs( const fvec &k, fvec *kdx, fvec *kdy, fvec *kdz )
{
	kdx->resize( k.size() );			kdy->resize( k.size() );				kdz->resize( k.size() );
	(*kdx)[0] =   k[2];					(*kdy)[0] =    k[3];					(*kdz)[0] =   k[1];
	(*kdx)[1] =   k[5];					(*kdy)[1] =    k[6];					(*kdz)[1] = 2*k[4];
	(*kdx)[2] = 2*k[7] - k[4];			(*kdy)[2] =    k[8];					(*kdz)[2] =   k[5];
	(*kdx)[3] =   k[8];					(*kdy)[3] = -2*k[7] - k[4];				(*kdz)[3] =   k[6];
	(*kdx)[4] =   k[10];				(*kdy)[4] =    k[11];					(*kdz)[4] = 3*k[9];
	(*kdx)[5] = 2*k[12] - 3*k[9];		(*kdy)[5] =    k[13];					(*kdz)[5] = 2*k[10];
	(*kdx)[6] =   k[13];				(*kdy)[6] = -2*k[12] - 3*k[9];			(*kdz)[6] = 2*k[11];
	(*kdx)[7] = 3*k[14] - 1/4*k[10];	(*kdy)[7] =  3*k[15] - 1/4*k[11];		(*kdz)[7] =   k[12];
	(*kdx)[8] = 6*k[15] - 1/2*k[11];	(*kdy)[8] = -6*k[14] - 1/2*k[10];		(*kdz)[8] =   k[13];
	(*kdx)[9] =   k[17];				(*kdy)[9] =    k[18];					(*kdz)[9] = 4*k[16];
	(*kdx)[10]= 2*k[19] - 6*k[16];		(*kdy)[10]=    k[20];					(*kdz)[10]= 3*k[17];
	(*kdx)[11]=   k[20];				(*kdy)[11]= -2*k[19] - 6*k[16];			(*kdz)[11]= 3*k[18];
	(*kdx)[12]= 3*k[21] - 3/4*k[17];	(*kdy)[12]=  3*k[22] - 3/4*k[18];		(*kdz)[12]= 2*k[19];
	(*kdx)[13]= 6*k[22] - 3/2*k[18];	(*kdy)[13]= -6*k[21] - 3/2*k[17];		(*kdz)[13]= 2*k[20];
	(*kdx)[14]= 4*k[23] - 1/6*k[19];	(*kdy)[13]=    k[24] + 1/12*k[20];		(*kdz)[14]=   k[21];
	(*kdx)[15]=   k[24] - 1/12*k[20];	(*kdy)[15]= -4*k[23] - 1/6*k[19];		(*kdz)[15]=   k[22];

	// Set the remaining coefficients to zero
	for ( unsigned i = 16; i != k.size(); ++i )
	{
		(*kdx)[i] = 0.0;
		(*kdy)[i] = 0.0;
		(*kdz)[i] = 0.0;
	}
}

/* 
 * FUNCTION: 	invertMatrix
 * DESRIPTION:	Calculates the inverse matrix of a 3x3 matrix.
 * INPUTS:		mat - matrix to be inverted.
 * OUTPUTS:		inv - inverse of the input matrix.
 */
float OptimalPositions::invertMatrix( const fmat3d &mat, fmat3d inv )
{
	float det = mat[0][0]*( mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1] ) -
				mat[0][1]*( mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0] ) +
				mat[0][2]*( mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0] );
	inv[0][0] = ( mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1] )/det;
	inv[0][1] = ( mat[0][2]*mat[2][1] - mat[0][1]*mat[2][2] )/det;
	inv[0][2] = ( mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1] )/det;
	inv[1][0] = ( mat[1][2]*mat[2][0] - mat[1][0]*mat[2][2] )/det;
	inv[1][1] = ( mat[0][0]*mat[2][2] - mat[0][2]*mat[2][0] )/det;
	inv[1][2] = ( mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2] )/det;
	inv[2][0] = ( mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0] )/det;
	inv[2][1] = ( mat[0][1]*mat[2][0] - mat[0][0]*mat[2][1] )/det;
	inv[2][2] = ( mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0] )/det;
	return det;
}
