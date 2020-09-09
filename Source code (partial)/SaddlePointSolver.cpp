#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_complex_math.h>

#include <cmath>
#include <iostream>
#include "SaddlePointSolver.hpp"
#include "RangeConverter.hpp"

///Global Period, Possibly better to define in main temporary way to change the range solver solves over
const double Pi = 3.1415926535897932384626433832795028841971693993751;
const double Period = 0.*Pi;

namespace
{

    struct MeritFunctionParams
    {
        MeritFunctionParams( const I_Action * f, double plong, double ptrans ) : myF( f ), myPLong( plong ), myPTrans( ptrans ) {};
        const I_Action * myF;
        double myPLong;
        double myPTrans;
    };


   int meritFunction_t1andt2 (const gsl_vector *x, void *params, gsl_vector* f)
    {

        MeritFunctionParams* p = static_cast< MeritFunctionParams* > ( params );

        const I_Action * myfunc = p->myF;

        /// get the derivative

        std::vector< gsl_complex > times;

        double t1real = RangeConverter::convertToBounded( gsl_vector_get( x, 0 ), Period , 3.*M_PI  + Period );
        double t2real = RangeConverter::convertToBounded( gsl_vector_get( x, 2 ), t1real+1e-11, 3.*M_PI + Period);
        double t1imag = gsl_vector_get( x, 1 );
        double t2imag = gsl_vector_get( x, 3 );

        times.push_back ( gsl_complex_rect( t1real , t1imag ) );
        times.push_back ( gsl_complex_rect( t2real , t2imag ) );
        times.push_back ( gsl_complex_rect( 0, 0 ) );

        gsl_complex deriv1 = myfunc->df( p->myPLong, p->myPTrans, p->myPLong, p->myPTrans, times, 0, false );
        gsl_complex deriv2 = myfunc->df( p->myPLong, p->myPTrans, p->myPLong, p->myPTrans, times, 1, false );



        /// trying to minimize both real and imaginary part of the function
        gsl_vector_set( f, 0 , GSL_REAL( deriv1 ) );
        gsl_vector_set( f, 1 , GSL_IMAG( deriv1 ) );
        gsl_vector_set( f, 2 , GSL_REAL( deriv2 ) );
        gsl_vector_set( f, 3 , GSL_IMAG( deriv2 ) );

        return GSL_SUCCESS;
    }


    int meritFunction_t3 (const gsl_vector *x, void *params, gsl_vector* f)
    {

        MeritFunctionParams* p = static_cast< MeritFunctionParams* > ( params );

        const I_Action * myfunc = p->myF;

        /// get the derivative

        std::vector< gsl_complex > times;

        double t3real = RangeConverter::convertToBounded( gsl_vector_get( x, 0 ), M_PI + Period, 3. * M_PI + Period);
        double t3imag = gsl_vector_get( x, 1 ) * gsl_vector_get( x, 1 );
        times.push_back( gsl_complex_rect( 0, 0 ) );
        times.push_back( gsl_complex_rect( 0, 0 ) );
        times.push_back ( gsl_complex_rect( t3real, t3imag ) );

        gsl_complex deriv3 = myfunc->df( p->myPLong, p->myPTrans, p->myPLong, p->myPTrans, times, 2, false);

        /// trying to minimize both real and imaginary part of the function


        gsl_vector_set( f, 0 , GSL_REAL( deriv3 ) );
        gsl_vector_set( f, 1 , GSL_IMAG( deriv3 ) );

        return GSL_SUCCESS;
    }



}


int SaddlePointSolver::solve( const I_Action * spf, double p_long, double p_trans, std::vector< gsl_complex >& times )
{

    const gsl_multiroot_fsolver_type * T = gsl_multiroot_fsolver_hybrid;
    gsl_multiroot_fsolver * s = gsl_multiroot_fsolver_alloc (T, 2 );

    gsl_multiroot_function f;
    MeritFunctionParams params ( spf, p_long, p_trans );


    f.f = &meritFunction_t3;
    f.params = &params;
    f.n = 2;

    gsl_vector* x = gsl_vector_alloc( 2 );

    double mappedrealval = RangeConverter::convertToUnbounded( GSL_REAL( times.at( 2 ) ), M_PI + Period, 3.*M_PI  + Period) ;
    gsl_vector_set( x, 0, mappedrealval );
    gsl_vector_set( x, 1, std::sqrt( GSL_IMAG( times.at( 2 ) ) ) );
    int status = gsl_multiroot_fsolver_set( s, &f, x );

    /// iterations

    int iter = 0;
    do{
        ++iter;
        status = gsl_multiroot_fsolver_iterate( s );
        ///std::cout << status << "\t" << gsl_vector_get( s->f, 0 ) << "\t" << gsl_vector_get( s->f, 1 ) << "\t" << gsl_vector_get( s->x, 0 ) << "\t" << gsl_vector_get( s->x, 1 ) << "\t" << std::endl;

        status = gsl_multiroot_test_residual( s->f, 1e-11 );
    }
    while( status == GSL_CONTINUE && iter < 1000 );

    if ( status == GSL_SUCCESS )
    {
            gsl_complex unmappedresult = gsl_complex_rect( RangeConverter::convertToBounded( gsl_vector_get( s->x, 0 ), M_PI + Period, 3.*M_PI + Period),
                                        gsl_vector_get( s->x, 1 )*gsl_vector_get( s->x, 1 ) );
            times.at( 2 ) = unmappedresult;

    }

    gsl_vector_free ( x );
    gsl_multiroot_fsolver_free( s );

    /// now solve for t1 and t2
    int status2  = 0;
    const gsl_multiroot_fsolver_type * T2 = gsl_multiroot_fsolver_hybrid;

    gsl_multiroot_fsolver * s2 = gsl_multiroot_fsolver_alloc ( T2, 4 );


    f.f = &meritFunction_t1andt2;
    f.params = &params;
    f.n = 4;
    gsl_vector* x2 = gsl_vector_alloc( 4 );

    mappedrealval = RangeConverter::convertToUnbounded( GSL_REAL ( times.at( 0 ) ), Period, 3.*M_PI  + Period);

    gsl_vector_set( x2, 0, mappedrealval );
    gsl_vector_set( x2, 1, GSL_IMAG( times.at( 0 ) ) );

    mappedrealval = RangeConverter::convertToUnbounded( GSL_REAL ( times.at( 1 ) ), GSL_REAL ( times.at( 0 ) ) + 1e-11, 3. * M_PI + Period);
    gsl_vector_set( x2, 2, mappedrealval );
    gsl_vector_set( x2, 3, GSL_IMAG( times.at( 1 ) ) );

    status2 = gsl_multiroot_fsolver_set( s2, &f, x2 );


     /// iterations

    iter = 0;
    do{
        ++iter;
        status2 = gsl_multiroot_fsolver_iterate( s2 );
        ///std::cout << status << "\t" << gsl_vector_get( s->f, 0 ) << "\t" << gsl_vector_get( s->f, 1 ) << "\t" << gsl_vector_get( s->x, 0 ) << "\t" << gsl_vector_get( s->x, 1 ) << "\t" << std::endl;


        status2 = gsl_multiroot_test_residual( s2->f, 1e-10 );
    }
    while( status2 == GSL_CONTINUE && iter < 1000 );

    if ( status2 == GSL_SUCCESS )
    {
            double t1real = RangeConverter::convertToBounded( gsl_vector_get( s2->x, 0 ), Period, 3. * M_PI + Period);
            gsl_complex unmappedresult = gsl_complex_rect( t1real, gsl_vector_get( s2->x, 1 ) );
            times.at( 0 ) = unmappedresult;

            double t2real = RangeConverter::convertToBounded( gsl_vector_get( s2->x, 2 ), t1real+1e-11, 3. * M_PI  + Period);

            unmappedresult = gsl_complex_rect( t2real, gsl_vector_get( s2->x, 3 ) );

            times.at( 1 ) = unmappedresult;
    }
    else
    {
    	//ToDo: Breaks here
    	std::cout << "Calibration 2 failed at iter=" << iter << std::endl;
        throw "whoopsie daisy";
    }

    gsl_vector_free ( x2 );
    gsl_multiroot_fsolver_free( s2 );

    return std::abs((double)status)+std::abs((double)status2);
}
