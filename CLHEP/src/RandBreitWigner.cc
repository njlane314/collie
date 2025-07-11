// $Id: RandBreitWigner.cc,v 1.1 2010/02/16 20:35:18 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                       --- RandBreitWigner ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Added methods to shoot arrays: 28th July 1997
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments: 16th Feb 1998
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/RandBreitWigner.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/DoubConv.hh"
#include <algorithm>	// for min() and max()
#include <cmath>

using namespace std;

namespace CLHEP {

std::string RandBreitWigner::name() const {return "RandBreitWigner";}
HepRandomEngine & RandBreitWigner::engine() {return *localEngine;}

RandBreitWigner::~RandBreitWigner() {
  if ( deleteEngine ) delete localEngine;
}

RandBreitWigner::RandBreitWigner(const RandBreitWigner& right)
 : defaultA(right.defaultA), defaultB(right.defaultB)
{;}

double RandBreitWigner::operator()() {
   return fire( defaultA, defaultB );
}

double RandBreitWigner::operator()( double a, double b ) {
   return fire( a, b );
}

double RandBreitWigner::operator()( double a, double b, double c ) {
   return fire( a, b, c );
}

double RandBreitWigner::shoot(double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*HepRandom::getTheEngine()->flat()-1.0;
   displ = 0.5*gamma*tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::shoot(double mean, double gamma, double cut)
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(2.0*cut/gamma);
   rval = 2.0*HepRandom::getTheEngine()->flat()-1.0;
   displ = 0.5*gamma*tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::shootM2(double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(-mean/gamma);
   rval = RandFlat::shoot(val, CLHEP::halfpi);
   displ = gamma*tan(rval);

   return sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::shootM2(double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = max(0.0,(mean-cut));
   lower = atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(lower, upper);
   displ = gamma*tan(rval);

   return sqrt(max(0.0, mean*mean + mean*displ));
}

void RandBreitWigner::shootArray ( const int size, double* vect )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( 1.0, 0.2 );
}

void RandBreitWigner::shootArray ( const int size, double* vect,
                                   double a, double b )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( a, b );
}

void RandBreitWigner::shootArray ( const int size, double* vect,
                                   double a, double b,
                                   double c )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( a, b, c );
}

//----------------

double RandBreitWigner::shoot(HepRandomEngine* anEngine,
                                 double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*anEngine->flat()-1.0;
   displ = 0.5*gamma*tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::shoot(HepRandomEngine* anEngine,
                                 double mean, double gamma, double cut )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(2.0*cut/gamma);
   rval = 2.0*anEngine->flat()-1.0;
   displ = 0.5*gamma*tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::shootM2(HepRandomEngine* anEngine,
                                   double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(-mean/gamma);
   rval = RandFlat::shoot(anEngine,val, CLHEP::halfpi);
   displ = gamma*tan(rval);

   return sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::shootM2(HepRandomEngine* anEngine,
                                   double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = max(0.0,(mean-cut));
   lower = atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(anEngine, lower, upper);
   displ = gamma*tan(rval);

   return sqrt( max(0.0, mean*mean+mean*displ) );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( anEngine, 1.0, 0.2 );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect,
                                   double a, double b )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( anEngine, a, b );
}

void RandBreitWigner::shootArray ( HepRandomEngine* anEngine,
                                   const int size, double* vect,
                                   double a, double b, double c )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot( anEngine, a, b, c );
}

//----------------

double RandBreitWigner::fire()
{
  return fire( defaultA, defaultB );
}

double RandBreitWigner::fire(double mean, double gamma)
{
   double rval, displ;

   rval = 2.0*localEngine->flat()-1.0;
   displ = 0.5*gamma*tan(rval*CLHEP::halfpi);

   return mean + displ;
}

double RandBreitWigner::fire(double mean, double gamma, double cut)
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(2.0*cut/gamma);
   rval = 2.0*localEngine->flat()-1.0;
   displ = 0.5*gamma*tan(rval*val);

   return mean + displ;
}

double RandBreitWigner::fireM2()
{
  return fireM2( defaultA, defaultB );
}

double RandBreitWigner::fireM2(double mean, double gamma )
{
   double val, rval, displ;

   if ( gamma == 0.0 ) return mean;
   val = atan(-mean/gamma);
   rval = RandFlat::shoot(localEngine,val, CLHEP::halfpi);
   displ = gamma*tan(rval);

   return sqrt(mean*mean + mean*displ);
}

double RandBreitWigner::fireM2(double mean, double gamma, double cut )
{
   double rval, displ;
   double lower, upper, tmp;

   if ( gamma == 0.0 ) return mean;
   tmp = max(0.0,(mean-cut));
   lower = atan( (tmp*tmp-mean*mean)/(mean*gamma) );
   upper = atan( ((mean+cut)*(mean+cut)-mean*mean)/(mean*gamma) );
   rval = RandFlat::shoot(localEngine,lower, upper);
   displ = gamma*tan(rval);

   return sqrt(max(0.0, mean*mean + mean*displ));
}

void RandBreitWigner::fireArray ( const int size, double* vect)
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire(defaultA, defaultB );
}

void RandBreitWigner::fireArray ( const int size, double* vect,
                                  double a, double b )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( a, b );
}

void RandBreitWigner::fireArray ( const int size, double* vect,
                                  double a, double b, double c )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( a, b, c );
}


std::ostream & RandBreitWigner::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultA);
  os << defaultA << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultB);
  os << defaultB << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
#ifdef REMOVED
  int pr=os.precision(20);
  os << " " << name() << "\n";
  os << defaultA << " " << defaultB << "\n";
  os.precision(pr);
  return os;
#endif
}

std::istream & RandBreitWigner::get ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != name()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read state of a "
    	      << name() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  if (possibleKeywordInput(is, "Uvec", defaultA)) {
    std::vector<unsigned long> t(2);
    is >> defaultA >> t[0] >> t[1]; defaultA = DoubConv::longs2double(t); 
    is >> defaultB >> t[0] >> t[1]; defaultB = DoubConv::longs2double(t); 
    return is;
  }
  // is >> defaultA encompassed by possibleKeywordInput
  is >> defaultB;
  return is;
}


}  // namespace CLHEP

