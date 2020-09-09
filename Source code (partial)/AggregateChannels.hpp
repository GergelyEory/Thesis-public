#ifndef __AGGREGATECHANNELS_HPP_
#define __AGGREGATECHANNELS_HPP_

#include <vector>
#include <fstream>

#include <gsl/gsl_complex.h>
#include "I_Action.hpp"




/// forward declaration
class ResultElement;

class AbstractPrefactor;

class AggregateChannels
{


public:
///typedef to shorten Amplitude notation
typedef std::vector<std::vector<std::vector <std::vector<gsl_complex> > > > AmpArray;
typedef std::vector <std::vector<gsl_complex> > AmpArray2;
typedef std::vector <std::vector<double> > ProbArray;

    ///constructors
AggregateChannels( const I_Action * f, const std::vector< std::vector< ResultElement > >& table_short
, const std::vector< std::vector< ResultElement > >& table_long, const std::vector<const AbstractPrefactor*> prefactors,std::vector<double> EnergyRangesIn);
    ///overload constructor so Amplitudes can be initialised
AggregateChannels( const I_Action * f, const std::vector< std::vector< ResultElement > >& table_short
, const std::vector< std::vector< ResultElement > >& table_long, const std::vector<const AbstractPrefactor*> prefactors, AmpArray Amplitude,std::vector<double> EnergyRangesIn);
    ///constructor so can load amplitudes values from file
AggregateChannels(int SizeLong, int SizeTrans, std::ostringstream& AmpFileName, const I_Action *f);
    ///constructor so can initialise amplitudes with specific values
AggregateChannels(int SizeLong, int SizeTrans, gsl_complex Initial, const I_Action *f,const std::vector< std::vector< ResultElement > >& table_short
, const std::vector< std::vector< ResultElement > >& table_long);


    ///Overloaded Operators
    AggregateChannels operator+(AggregateChannels const &other);
    AggregateChannels operator*(const gsl_complex other) const;


    ///Public Functions
    int Mul(gsl_complex Coef);

    int LoadTimes() ;
    int LoadAmplitudesFromFile(std::ostringstream& AmpFileName);
    int SaveAmplitudes(bool PrintCount, std::ostream& oFile, bool printToFile);
    int SetAmplitudes(gsl_complex);

    int IntegratePerp();
    int IntegratePerp(std::ostringstream &OutFile);
    int FixPerp(double Value, double max_energy, double increment);
    int FixPerp(std::ostringstream &OutFile, double Value, double max_energy, double increment);
    int IntegratePar();
    int IntegratePar(std::ostringstream &OutFile);
    double ProbabilityMean();
    int extractP1Partial();
    int extractP2Partial();
    int extractP1Partial(std::ostringstream &FileNameOut);
    int extractP2Partial(std::ostringstream &FileNameOut);
    int extractP1Angle();
    int extractP2Angle();
    int extractPparAngle();
    int extractP1Angle(std::ostringstream &FileNameOut);
    int extractP2Angle(std::ostringstream &FileNameOut);
    int extractPparAngle(std::ostringstream &FileNameOut);
    int IncoherentSymmetrise(std::vector<bool> MomentumSwap);
    int CoherentSymmetrise(std::vector<bool> MomentumSwap);
    int IncoherentSymmetrise(std::vector<bool> MomentumSwap, std::ostringstream &OutFile);

    void Add2_Prob(const ProbArray &P1);


    const AmpArray &getAmplitude() const;
    const ProbArray &getProbability() const;
    const AmpArray2 &getPartialP1Amplitude() const;
    const ProbArray &getPartialP1Probability() const;
    const AmpArray2 &getPartialP2Amplitude() const;
    const ProbArray &getPartialP2Probability() const;

    int clearProbability();
    int clearSymmetrisedProbability();
    int clearAmplitude();
    int clearP1Partial();
    int clearP2Partial();
    int clearP1Angle();
    int clearP2Angle();
    int clearPparAngle();


    int SetProbability(ProbArray &Prob);

    int PrintArray2D(const ProbArray& Array2D, std::ostream& oFile, bool GnuPlot) const;

    int PrintAmplitudeToFile(std::ostringstream& AmpFileName) const;
    int PrintProbabilityToFile(std::ostream& longstream, bool ExtraLine) const;
    int PrintAmpArray2D(const AmpArray2& Array2D, std::ostream& oFile, bool GnuPlot) const;
    int PrintP1PartialToFile(std::ostream& longstream, bool GnuPlot) const;
    int PrintP2PartialToFile(std::ostream& longstream, bool GnuPlot) const;
    int PrintSymmetrisedProbability(std::ostream& longstream, bool ExtraLine) const;


private:

    ///Private Functions
    AmpArray Amplitude_Add(const AmpArray &A1, const AmpArray &A2) const;
    AmpArray Scalar_Mul(const AmpArray &A1, const gsl_complex Coef) const;

    /// getPrefactor12
    void getPrefactor12( std::vector< gsl_complex >& retval, const double p1long, const double p1trans, const std::vector< gsl_complex >& times ) const;
    /// getPrefactor3
    gsl_complex getPrefactor3( const std::vector<gsl_complex>& prefactor12result, const double p1long, const double p1trans, const double p2long, const double p2trans, const std::vector< gsl_complex >& times ) const;

    ///Private Data
    const I_Action * myF;
    std::vector< std::vector< ResultElement > > myTableShort;
    std::vector< std::vector< ResultElement > > myTableLong;
    std::vector< const AbstractPrefactor * > myPrefactors;

    AmpArray TransitionAmplitude;
    AmpArray2 TransitionAmplitudeE1;

    std::vector<double > EnergyRanges;
    ProbArray TransitionProbability;

    AmpArray2 P1PartialAmplitude;
    AmpArray2 P2PartialAmplitude;
    ProbArray P1PartialProbability;
    ProbArray P2PartialProbability;
    ProbArray P1Angle;
    ProbArray P2Angle;
    ProbArray PparAngle;




    ProbArray SymmetrisedTransitionProbability;
    std::vector< std::vector< ResultElement > > cached_t1t2;
    std::vector< std::vector< ResultElement > > cached_t3;
    AmpArray2 d33_matrix;
    AmpArray2 actionvalue_t3_matrix;


};

///NonClass operators
AggregateChannels operator*(gsl_complex other, AggregateChannels Aggregator);

#endif // AGGREGATECHANNELS_HPP_INCLUDED
