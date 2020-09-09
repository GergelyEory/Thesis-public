#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <iomanip>
#include <exception>

#include "AggregateChannels.hpp"
#include "ActionAnalytics.hpp"
#include "I_Action.hpp"
#include "ResultGrid.hpp"
#include "AbstractPrefactor.hpp"
#include "Xtra_Funcs.hpp"

const double Pi = 3.1415926535897932384626433832795028841971693993751;
const gsl_complex I = gsl_complex_rect(0,1);

typedef std::vector<std::vector<std::vector<std::vector<gsl_complex> > > > AmpArray;
typedef std::vector<std::vector<std::vector<gsl_complex> > > AmpArrayRow;
typedef std::vector<std::vector<gsl_complex> >AmpArray2;
typedef std::vector<gsl_complex>AmpArray2Row;
typedef std::vector <std::vector<double> > ProbArray;
typedef std::vector<double> ProbArrayRow;
typedef std::vector <std::vector<ResultElement> > ResultArray;
typedef std::vector <ResultElement>ResultArrayRow;

///Constructors for class

///Contruct Aggregator Ready to Calculate Amplitudes
AggregateChannels::AggregateChannels(  const I_Action * f, const std::vector< std::vector< ResultElement > >& table_short
, const std::vector< std::vector< ResultElement > >& table_long, const std::vector<const AbstractPrefactor*> prefactors, std::vector<double> EnergyRangesIn)
: myF( f ), myTableShort( table_short ), myTableLong( table_long ), myPrefactors (prefactors ), EnergyRanges(EnergyRangesIn)
{}

///Contruct Aggregator with Amplitudes Already Set.
AggregateChannels::AggregateChannels(  const I_Action * f, const std::vector< std::vector< ResultElement > >& table_short
, const std::vector< std::vector< ResultElement > >& table_long, const std::vector<const AbstractPrefactor*> prefactors, AmpArray Amplitude, std::vector<double> EnergyRangesIn)
: myF( f ), myTableShort( table_short ), myTableLong( table_long ), myPrefactors (prefactors ), TransitionAmplitude(Amplitude), EnergyRanges(EnergyRangesIn)
{
    LoadTimes();
}


///Contruct Aggregator with Amplitudes set to a value and no other variables set
AggregateChannels::AggregateChannels(int SizeLong, int SizeTrans, gsl_complex Initial, const I_Action *f
,const std::vector< std::vector< ResultElement > >& table_short, const std::vector< std::vector< ResultElement > >& table_long): myF(f), myTableShort(table_short), myTableLong(table_long)
{
    LoadTimes();
    SizeTrans=std::round(SizeTrans/2.);
    std::vector<gsl_complex> Times_Zero(3,gsl_complex_rect(0,0));
    ResultElement ZeroEl(Times_Zero, 0, 0.0, 0.0);
    myTableShort.resize(SizeLong,std::vector< ResultElement >(SizeTrans, ZeroEl));
    myTableLong.resize(SizeLong,std::vector< ResultElement >(SizeTrans, ZeroEl));

    TransitionAmplitude.assign(SizeLong,AmpArrayRow(SizeLong,AmpArray2(SizeTrans,AmpArray2Row(SizeTrans,Initial))));
}


///Operator Overloading
AggregateChannels AggregateChannels::operator+(AggregateChannels const &other)
{
    ///Could have checks here to prevent class abuse

    AmpArray A2 = other.getAmplitude();

    return AggregateChannels(myF, myTableShort, myTableLong, myPrefactors, Amplitude_Add(TransitionAmplitude, A2), EnergyRanges);
}

AggregateChannels AggregateChannels::operator*(const gsl_complex other) const
{
    AmpArray NewAmplitude = Scalar_Mul(TransitionAmplitude,other);
    return AggregateChannels(myF, myTableShort, myTableLong, myPrefactors, NewAmplitude, EnergyRanges);
}
///Commutativity for scalar multiplication
AggregateChannels operator*(gsl_complex other, AggregateChannels Aggregator)
{
    return Aggregator*other;
}

int AggregateChannels::Mul(gsl_complex Coef)
{
    TransitionAmplitude = Scalar_Mul(TransitionAmplitude, Coef);
    return 0;
}

AmpArray AggregateChannels::Amplitude_Add(const AmpArray &A1, const AmpArray &A2) const
{
    AmpArray Aout;

    AmpArray::const_iterator it1=A1.begin();
    AmpArray::const_iterator it2=A2.begin();

    if((A1.size()!=A2.size())||(A1.at(0).size()!=A2.at(0).size())||(A1.at(0).at(0).size()!=A2.at(0).at(0).size())||(A1.at(0).at(0).at(0).size()!=A2.at(0).at(0).at(0).size()))
    {std::cerr<<"Error Amplitudes Added are Different Sizes!\n A Zero Sized Amplitude will be returned\n";
    std::cout<<A1.size()<<", "<<A2.size()<<std::endl;
    std::cout<<A1.at(0).size()<<", "<<A2.at(0).size()<<std::endl;
    std::cout<<A1.at(0).at(0).size()<<", "<<A2.at(0).at(0).size()<<std::endl;
    std::cout<<A1.at(0).at(0).at(0).size()<<", "<<A2.at(0).at(0).at(0).size()<<std::endl;
    return Aout;}

    for(int i=0;it2!=A2.end(); ++it1, ++it2,++i)
    {
        AmpArrayRow Row;
        AmpArrayRow::const_iterator it12=it1->begin();
        AmpArrayRow::const_iterator it22=it2->begin();
        for(int j=0;it22!=it2->end(); ++it12, ++it22,++j)
        {
            AmpArray2 Row2;
            AmpArray2::const_iterator it13=it12->begin();
            AmpArray2::const_iterator it23=it22->begin();
            for(int ii=0;it23!=it22->end(); ++it13, ++it23,++ii)
            {
                AmpArray2Row Row3;
                AmpArray2Row::const_iterator it14=it13->begin();
                AmpArray2Row::const_iterator it24=it23->begin();

                for(int jj=0;it24!=it23->end(); ++it14, ++it24, ++jj)
                {

                    gsl_complex tmp = gsl_complex_add(*it14, *it24);
                    Row3.push_back(tmp);
                }
                Row2.push_back(Row3);
                Row3.clear();
            }
            Row.push_back(Row2);
            Row2.clear();
        }
        Aout.push_back(Row);
        Row.clear();
    }
    return Aout;

}

AmpArray AggregateChannels::Scalar_Mul(const AmpArray &A1, const gsl_complex Coef) const
{
    AmpArray NewAmplitude;
    if(GSL_REAL(Coef)==1&&GSL_IMAG(Coef)==0)
    {
        NewAmplitude = A1;
    }
    else
    {
        AmpArray::const_iterator it=A1.begin();
        for(int i=0;it!=A1.end(); ++it,++i)
        {
            AmpArrayRow NewAmpRow;
            AmpArrayRow::const_iterator it2 = it->begin();
            for(int j=0;it2!=it->end(); ++it2,++j)
            {
                AmpArray2 NewAmp2;
                AmpArray2::const_iterator it3=it2->begin();
                for(int ii=0;it3!=it2->end(); ++it3, ++ii)
                {
                    std::vector<gsl_complex> NewAmp2Row;
                    std::vector<gsl_complex>::const_iterator it4 = it3->begin();
                    for(int jj=0;it4!=it3->end(); ++it4, ++jj)
                    {
                        gsl_complex tmp = gsl_complex_mul(*it4,Coef);
                        NewAmp2Row.push_back(tmp);
                    }
                    NewAmp2.push_back(NewAmp2Row);
                    NewAmp2Row.clear();
                }
                NewAmpRow.push_back(NewAmp2);
                NewAmp2.clear();
            }
            NewAmplitude.push_back(NewAmpRow);
            NewAmpRow.clear();
        }
        //std::cout<<GSL2C(A1.at(10).at(10).at(10).at(10))<<" -> "<<GSL2C(NewAmplitude.at(10).at(10).at(10).at(10))<<std::endl;
    }

    return NewAmplitude;
}

void AggregateChannels::Add2_Prob(const ProbArray &P1)
{

    ProbArray::iterator it1 = TransitionProbability.begin();
    ProbArray::const_iterator itP1 = P1.begin();

    for(int i=0 ;it1!=TransitionProbability.end(); ++it1 , ++itP1, ++i)
    {
        std::vector<double>::iterator it2 = it1->begin();
        std::vector<double>::const_iterator itP2 = itP1->begin();
        for(int j=0 ;it2!=it1->end(); ++it2, ++itP2, ++j)
        {
            ///std::cout<<i<<", "<<j<<std::endl;
            *it2 = *it2 + *itP2;

        }
    }
}

///Loads previously calculated times into class. Fills cached_t1t2 and cached_t3 has capability to invert/ flip times
int AggregateChannels::LoadTimes()
{
        /// loop over the files (taba, tabb, tabc, tabd)
        ActionAnalytics a( myF );
        std::vector< double > row_vector;
        //std::cout << "caching t3 table...\n";
        std::cout.flush();

        std::vector< ResultElement > tempvect1t2;

        // Iterators over radial spokes from first to last
        std::vector< std::vector< ResultElement > >::const_iterator short_spokes_it = myTableShort.begin();
        std::vector< std::vector< ResultElement > >::const_iterator long_spokes_it = myTableLong.begin();


        if(myTableShort.size()!=myTableLong.size()){std::cerr<<"ErrorShort and Long not same size"<<std::endl;return -1;}
        for (int i=0 ; short_spokes_it != myTableShort.end() ; ++short_spokes_it, ++long_spokes_it, ++i)
        {
        	// Iterators over elements in each radial spoke from first to last
            std::vector< ResultElement >::const_iterator short_radial_it = short_spokes_it->begin();
            std::vector< ResultElement >::const_iterator long_radial_it = long_spokes_it->begin();

            for (int j=0 ; short_radial_it != short_spokes_it->end() ; ++short_radial_it, ++long_radial_it, ++j )
            {

                double p2long = short_radial_it->getPLong();
                double p2trans = short_radial_it->getPTrans();

                std::vector< gsl_complex > valuest1t2;

				valuest1t2.push_back( short_radial_it->getValues()[ 0 ] );
				valuest1t2.push_back( long_radial_it->getValues()[ 0 ] );

				valuest1t2.push_back( short_radial_it->getValues()[ 1 ] );
				valuest1t2.push_back( long_radial_it->getValues()[ 1 ] );



                ResultElement elemt1t2( valuest1t2, false, p2long, p2trans );
                tempvect1t2.push_back( elemt1t2 );
            }

            cached_t1t2.push_back ( tempvect1t2 );
            tempvect1t2.clear();
        }

        return 0;

}

///Generates a npar X nperp X npar X nperp of all the complex amplitudes. Prefactor are included if specified in here.
///HalfPerp is called as true for some bloody reason
///Set HalfPerp to False as radial grid might not be exactly symmetric
int AggregateChannels::SaveAmplitudes(bool PrintCount, std::ostream& oFile, bool printToFile) {

    TransitionAmplitudeE1.reserve(cached_t1t2.size());

    // Iterator over grid of ResultElement thetas, their values contain t1 and t2 data
    ResultArray::const_iterator it = cached_t1t2.begin();
    double percent=0;
    for(int i=0; it!=cached_t1t2.end(); ++it, ++i) {
    	//Percentage Counter
        if(PrintCount) {

        	percent = 100*(double)(i+1) / myTableShort.capacity();
			std::cerr<<"\r"<<"Completed "<<std::setprecision(3)<<percent<<"% of Calculating Amplitudes";
			std::cerr.flush();
        }

        ActionAnalytics a( myF );
		std::vector< ResultElement >::const_iterator it3 = it->begin();

		std::vector<gsl_complex> TransitionAmplitudeSpoke;
		TransitionAmplitudeSpoke.reserve(it->size());

		// Iterating over the radial elements
		for (int ii = 0; it3 != it->end(); ++it3, ++ii) {
			double p1long = it3->getPLong();
			double p1trans = it3->getPTrans();
			std::vector< gsl_complex > times_t1t2 = it3->getValues();

			std::vector< gsl_complex > times_short;
			std::vector< gsl_complex > times_long;

			/// short times first
			times_short.push_back( times_t1t2[ 0 ] );
			times_short.push_back( times_t1t2[ 2 ] );
			times_short.push_back( gsl_complex_rect( 0,0 ) );

			/// the the long times
			times_long.push_back( times_t1t2[ 1 ] );
			times_long.push_back( times_t1t2[ 3 ] );
			times_long.push_back( gsl_complex_rect( 0,0 ) );

			gsl_complex action_value12_long = myF->f12( p1long, p1trans, times_long[ 0 ], times_long[ 1 ] );
			gsl_complex action_value12_short = myF->f12( p1long, p1trans, times_short[ 0 ], times_short[ 1 ] );

			// Check for Stoke transition
			bool stoke12 = a.StokesTransition12(p1long, p1trans, times_short, times_long);

			// Matrix element, without prefactors
			gsl_complex Mij;

			//Need to define but not calculated, set to 1.+0.i
			gsl_complex Along = gsl_complex_rect(1,0);
			gsl_complex Ashort = gsl_complex_rect(1,0);

			// Changed to use the 12 action, t3 related action and any prefactors are not calculated.
			if (stoke12)
			{
				Mij =  a.getMij_C(action_value12_short, action_value12_long, Ashort, Along);
			}
			else
			{
				Mij =  a.getMij_Q(action_value12_long, action_value12_short, Along, Ashort);
			}

			if(printToFile) {
				//std::cout<< "Printing amplitudes to file" << std::endl;
				oFile<<std::setprecision(12) <<  p1long   << " " << p1trans  <<" "<< GSL_REAL(Mij) << " " << GSL_IMAG(Mij) << std::endl;
			}

			TransitionAmplitudeSpoke.push_back(Mij);
		}
		TransitionAmplitudeE1.push_back(TransitionAmplitudeSpoke);

    }
    std::cerr<<std::endl;
    if(TransitionAmplitudeE1.size()!=0) {
    	return 0;
    } else {
    	return -1;
    }

}

///Function to return private class variables
const AmpArray &AggregateChannels::getAmplitude() const
{
    return TransitionAmplitude;
}
const ProbArray &AggregateChannels::getProbability() const
{
    return TransitionProbability;
}

///Set TransionProbabnilty with specific array
int AggregateChannels::SetProbability(ProbArray &Prob)
{
    TransitionProbability = Prob;
    return 0;
}


///Print Data to human-readable text files
int AggregateChannels::PrintArray2D(const ProbArray& Array2D, std::ostream& oFile, bool GnuPlot) const
{
    ProbArray::const_iterator it1 = Array2D.begin();
    for(int i=0; it1!=Array2D.end(); ++it1, ++i)
    {
        std::vector<double >::const_iterator it2 = it1->begin();
            for (int j=0;it2 != it1->end(); ++it2, ++j)
            {
                oFile<<std::setprecision(12) <<  i   << " " << j  <<" "<< *it2 << std::endl;
            }

         if(GnuPlot==true){oFile << std::endl;}
    }

    return 0;
}


int AggregateChannels::PrintProbabilityToFile(std::ostream& longstream, bool ExtraLine) const
{
    PrintArray2D(TransitionProbability, longstream, ExtraLine);
    return 0;
}

int AggregateChannels::PrintSymmetrisedProbability(std::ostream& longstream, bool ExtraLine) const
{
    PrintArray2D(SymmetrisedTransitionProbability, longstream, ExtraLine);
    return 0;
}

int AggregateChannels::PrintP1PartialToFile(std::ostream& longstream, bool GnuPlot) const
{
    PrintArray2D(P1PartialProbability, longstream, GnuPlot);
    return 0;
}
int AggregateChannels::PrintP2PartialToFile(std::ostream& longstream, bool GnuPlot) const
{
    PrintArray2D(P2PartialProbability, longstream, GnuPlot);
    return 0;
}




