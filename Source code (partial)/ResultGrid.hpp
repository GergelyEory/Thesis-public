#ifndef __RESULTGRID_HPP__
#define __RESULTGRID_HPP__

#include "ResultElement.hpp"
#include "I_Action.hpp"

#include <vector>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>




class ResultGrid
{
public:
    /// ctor
    ResultGrid (    const I_Action * f   			/// saddlepointfunction
                ,   double rStartIn
                ,   double rEndIn
                ,   double thetaStartIn
                ,   double thetaEndIn
				,	double rStepIn
				,	double thetaStepIn
				,	int Nrad
				,	int Ntheta
                ,   const std::vector< gsl_complex >& startguess /// start guess
                );
    /// populate table
    int populateTable();
    int PiTransformTable();
    int ReverseMomenta();

    ///dtor
    ~ResultGrid()
    {
        /// DO NOT DELETE THE POINTER - it's owned by the parent class !!!
        ///delete myAction;
    }

    /// solve
    ResultElement solvePoint(int iTheta, int iRad, int direction) const;


    /// get the results
    const std::vector< std::vector< ResultElement > >& getResultTable() const{ return myResultTable; }

private:
    const I_Action * myAction;
    double rStart;
    double rEnd;
    double thetaStart;
    double thetaEnd;
    double rStep;
    double thetaStep;
    int Nrad;
    int Ntheta;
    std::vector< gsl_complex > myStartGuess;
    std::vector< std::vector< ResultElement > > myResultTable;
};


#endif



