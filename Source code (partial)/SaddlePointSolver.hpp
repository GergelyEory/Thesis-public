/* saddlepoint solver class */

#ifndef __SADDLEPOINTSOLVER_HPP__
#define __SADDLEPOINTSOLVER_HPP__

#include <vector>
#include <gsl/gsl_complex.h>
#include "I_Action.hpp"


class SaddlePointSolver
{
public:
    /// ctor
    SaddlePointSolver( const I_Action * f ) : myF ( f ){};

    /// solve for some t - note times passed by non-const ref because some of them can (and will) be modified
    int solve( const I_Action * spf, double p_long, double p_trans, std::vector< gsl_complex >& times );

private:
    /// the function to be solved
    const I_Action * myF;


};


#endif



