#ifndef __RESULTELEMENT_HPP__
#define __RESULTELEMENT_HPP__

#include <gsl/gsl_complex.h>

#include <vector>

class ResultElement
{
public:
    /// dummy ctor
    ResultElement(){}; // @suppress("Class members should be properly initialized")
    /// ctor
    ResultElement ( const std::vector< gsl_complex >& values,
                    bool interpolated,
                    double plong,
                    double ptrans )
                    : myValues ( values )
                    , myInterpolated ( interpolated )
                    , myPLong( plong )
                    , myPTrans( ptrans )
                    {};
    const std::vector< gsl_complex >& getValues() const
    {
            return myValues;
    }

    void setValues(std::vector<gsl_complex>& Values);

    bool isInterpolated() const
    {
            return myInterpolated;
    }

    double getPLong() const
    {
            return myPLong;
    }

    double getPTrans() const
    {
            return myPTrans;
    }

    void addPi();

    void ReverseMomenta();



private:
   std::vector< gsl_complex > myValues;
   bool myInterpolated;
   double myPLong;
   double myPTrans;

 };



#endif /// end header guard
