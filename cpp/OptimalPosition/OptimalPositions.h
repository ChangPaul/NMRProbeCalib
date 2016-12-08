/*
 * OptimalPositions.h
 *
 *  Created on: Jun 11, 2015
 *      Author: Paul Chang
 */

#ifndef OPTIMALPOSITIONS_H_
#define OPTIMALPOSITIONS_H_

#include "../spha.h"
#include <vector>
#include <iostream>

using namespace std;

//==================================================================
//                        Type Definitions
//==================================================================

#define ALGORITHM 0			// 0: Newton's, 1: Gradient descent

typedef float fvec3d[3];
typedef float fmat3d[3][3];
typedef vector< float > fvec;
#define abs(x) ((x) < 0 ? -(x) : (x))

//==================================================================
//                       Class Declaration
//==================================================================

class OptimalPositions
{

public:
	// Constructor and Destructor
	OptimalPositions( const fvec &xgrad, const fvec &ygrad, const fvec &zgrad );
	virtual ~OptimalPositions() {}
	
	// Calculation method
	float Calculate( fvec3d freq );

private:
	float OptimStep		( const fvec3d &pos, const fvec3d &freq, fvec3d step );
	void  EstimateError	( const fvec3d &pos, const fvec3d &freq, fvec3d err );
	inline float Fval	( const fvec3d &pos, const fvec3d &freq );
	void  dFval			( const fvec3d &pos, const fvec3d &freq, fvec3d dF );
	void  ddFval		( const fvec3d &pos, const fvec3d &freq, fmat3d ddF );
	void  dFdx			( const fvec3d &pos, const vector< float > &k, fvec3d dx );
	void  derivCoeffs	( const fvec &k, fvec *kdx, fvec *kdy, fvec *kdz );
	float invertMatrix	( const fmat3d &mat, fmat3d inv );

private:
	int num_fn;
	vector< float > kx, ky, kz;
};

#endif /* OPTIMALPOSITIONS_H_ */
