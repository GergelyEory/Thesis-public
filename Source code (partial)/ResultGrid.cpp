#include "ResultGrid.hpp"
#include "SaddlePointSolver.hpp"

#include <cmath>
#include <iostream>
#include <algorithm>

//ctor
ResultGrid::ResultGrid(		const I_Action *f,
							double rStartIn,
							double rEndIn,
							double thetaStartIn,
							double thetaEndIn,
							double rStepIn,
							double thetaStepIn,
							int NradIn,
							int NthetaIn,
							const std::vector<gsl_complex>& startguess
							) :
							myAction (f),
							rStart (rStartIn),
							rEnd (rEndIn),
							thetaStart (thetaStartIn),
							thetaEnd (thetaEndIn),
							rStep (rStepIn),
							thetaStep (thetaStepIn),
							Nrad (NradIn),
							Ntheta (NthetaIn),
							myStartGuess (startguess)
{}

bool testing = false;

int ResultGrid::populateTable() {

	int iTheta = 0;

	while (iTheta < Ntheta) {
		myResultTable.push_back(std::vector<ResultElement> () );
		std::vector<ResultElement>& currentSpoke = myResultTable.at(iTheta);

		// Direction along the spoke. 1 is out, -1 is in
		int directionSpoke = ((iTheta % 2) == 0 ? 1 : -1);

		// Radial index starting point
		int iRad = (directionSpoke > 0 ? 0 : Nrad -1);

		while ((iRad >= 0) && (iRad < Nrad)) {
			if(testing) {
				std::cout<<"trying to solve elem"<<std::endl;
				std::cout<<"iRad = "<<iRad<<" iTheta = "<<iTheta<<std::endl;
			}
			ResultElement elem = solvePoint(iTheta, iRad, directionSpoke);
			if(testing) {
				std::cout<<"elem solved"<<std::endl;
			}
			// Failed to solve for the point, exit.
			if (elem.getValues().size() == 0) {
				throw elem;
				// Failed to fill table. Exit code -1
				return -1;
			}

			if(directionSpoke > 0) {
				currentSpoke.push_back(elem);
			} else {
				std::vector< ResultElement >::iterator elem_it = currentSpoke.begin();
				currentSpoke.insert( elem_it, elem );
			}

			// Set start guess to last solved point
			myStartGuess = elem.getValues();

			iRad += directionSpoke;
		}
		iTheta++;
	}
	// Successfully filled table. Exit code 0
	return 0;
}

ResultElement addPiEl(ResultElement tempEl) {
	tempEl.addPi();
	return tempEl;
}

std::vector<ResultElement> addPiVec(std::vector<ResultElement>ResultVec) {
   std::transform(ResultVec.begin(), ResultVec.end(), ResultVec.begin(), addPiEl);
   return ResultVec;
}

int ResultGrid::PiTransformTable() {
	std::transform(myResultTable.begin(), myResultTable.end(), myResultTable.begin(), addPiVec);
	return 0;
}

ResultElement ReverseMomentaEl(ResultElement TempEl) {
    TempEl.ReverseMomenta();
    return TempEl;
}

std::vector<ResultElement> ReverseMomentaVec( std::vector<ResultElement> ResultVec) {
    std::transform(ResultVec.begin(), ResultVec.end(),ResultVec.begin(), ReverseMomentaEl);
    return ResultVec;
}

int ResultGrid::ReverseMomenta() {
	std::transform(myResultTable.begin(), myResultTable.end(), myResultTable.begin(), ReverseMomentaVec);
	return 0;
}

ResultElement ResultGrid::solvePoint(int iTheta, int iRad, int direction) const {
	SaddlePointSolver solver(myAction);
	// Check point is a nearby previously solved value
	int checkTheta = ( iTheta == 0 ? 0 : iTheta -1);

	int checkRad = iRad;
	if ((direction > 0) && iRad > 0) {
		checkRad = iRad -1;
	} else if ((direction < 0) && (iRad < Nrad -1)) {
		checkRad = iRad + 1;
	}

	double prevTheta = thetaStart + (thetaStep * checkTheta);
	double prevRad = rStart + (rStep * checkRad);
	double prevLong = prevRad * std::cos(prevTheta);
	double prevTrans = prevRad * std::sin(prevTheta);

	double currTheta = thetaStart + (thetaStep * iTheta);
	double currRad = rStart + (rStep * iRad);
	double currLong = currRad * std::cos(currTheta);
	double currTrans = currRad * std::sin(currTheta);

	double tempLong = currLong;
	double tempTrans = currTrans;


	/// set the guessvalue if we are not in the starting point

	std::vector< gsl_complex > startguess = myStartGuess;
	if (iTheta + iRad > 0) {
		startguess = myResultTable.at(checkTheta).at(checkRad).getValues();
	}

	int maxDepth = 100;
	int currentDepth = 0;
	bool reducedStep = false;
	bool statusSuccess = false;


	while (!statusSuccess && (currentDepth < maxDepth)) {

		if (currentDepth == 0 && (iTheta + iRad > 0)) {
			// Reset startguess to check point when at depth 0
			startguess = myResultTable.at(checkTheta).at(checkRad).getValues();
		}

		std::vector<gsl_complex> estimate = startguess;
		//solve point (iTheta, iRad)
		int statusOut = solver.solve(myAction, tempLong, tempTrans, estimate);

		// Solver changes the time values, so need to create new variables? not quite sure but keeping with it as i dont wanna break it
		std::vector<gsl_complex> returnvalue = estimate;
		//solve check point (checkTheta, checkRad)
		int statusReturn = solver.solve(myAction, prevLong, prevTrans, returnvalue);

		double distancefromguess = 0;
		if (iTheta + iRad > 1) {
			for ( int iii = 0 ; iii < 3 ; ++iii ) {
				distancefromguess += gsl_complex_abs ( gsl_complex_sub( startguess.at( iii ), returnvalue.at( iii ) ) );
			}
		}

		// If one of the solvers failed or distance from guess is too large
		if ((statusOut != 0 || statusReturn != 0) || (iTheta + iRad > 0 && distancefromguess > 1e-6)) {
			// bisect step between points, attempting to solve smaller step
			tempLong = (tempLong + prevLong) /2;
			tempTrans = (tempTrans + prevTrans) / 2;
			reducedStep = true;
			estimate = startguess;
		} else {
			//we have a match
			if (checkTheta != 0) {
				startguess = estimate;
			}
			if (!reducedStep) {
				//break bisection
				statusSuccess = true;
				return ResultElement(estimate, false, currLong, currTrans);
			} else {
				// try again without reduced step with new estimate
				reducedStep = false;
				prevLong = tempLong;
				prevTrans = tempTrans;
				// reset temp variable to current point
				tempLong = currLong;
				tempTrans = currTrans;
				startguess = estimate;
			}
		}
		++currentDepth;

		if (currentDepth > 5) {
			std::cout<<GSL_REAL(startguess.at(2))<<", "<<GSL_IMAG(startguess.at(2))<<std::endl;
			std::cout<<"i = "<<iTheta+1<<", j = "<<iRad+1<<"\t"<<prevLong<<", "<<prevTrans<<"\t"<<distancefromguess<<std::endl;
		}

		// Return dummy ResultElement if solver unsuccessful
	}

	return ResultElement();

}
