#include "ResultElement.hpp"

#include<vector>
#include<algorithm>
#include <gsl/gsl_complex_math.h>
const double Pi = 3.1415926535897932384626433832795028841971693993751;

gsl_complex gsl_complex_add_pi( gsl_complex z){return gsl_complex_add_real(z,Pi);}

void ResultElement::addPi()
{
    std::transform(myValues.begin(), myValues.end(), myValues.begin(), gsl_complex_add_pi);
}

void ResultElement::setValues(std::vector<gsl_complex>& Values)
{
    myValues = Values;
}

void ResultElement::ReverseMomenta()
{
    myPLong=-myPLong;
    myPTrans=-myPTrans;
}
