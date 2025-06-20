// $Id: RandGauss.cc,v 1.1 2010/02/16 20:35:18 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandGauss ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Added methods to shoot arrays: 28th July 1997
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments. Introduced method normal()
//                  for computation in fire(): 16th Feb 1998
// Gabriele Cosmo - Relocated static data from HepRandom: 5th Jan 1999
// M Fischler     - Copy constructor should supply right engine to HepRandom:
//		    1/26/00.
// M Fischler     - Workaround for problem of non-reproducing saveEngineStatus
//		    by saving cached gaussian.  March 2000.
// M Fischler     - Avoiding hang when file not found in restoreEngineStatus 
//                  12/3/04
// M Fischler     - put and get to/from streams 12/8/04
// M Fischler     - save and restore dist to streams 12/20/04
// M Fischler	  - put/get to/from streams uses pairs of ulongs when
//		    storing doubles avoid problems with precision.
//		    Similarly for saveEngineStatus and RestoreEngineStatus
//		    and for save/restore distState
//		    Care was taken that old-form output can still be read back.
//			4/14/05
//              
// =======================================================================

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/DoubConv.hh"
#include <string.h>
#include <cmath>	// for log()

namespace CLHEP {

std::string RandGauss::name() const {return "RandGauss";}
HepRandomEngine & RandGauss::engine() {return *localEngine;}

// Initialisation of static data
bool RandGauss::set_st = false;
double RandGauss::nextGauss_st = 0.0;

RandGauss::~RandGauss() {
  if ( deleteEngine ) delete localEngine;
}

RandGauss::RandGauss(const RandGauss& right)
 : HepRandom(right.getTheEngine()), 
   defaultMean(right.defaultMean), defaultStdDev(right.defaultStdDev)
{;}

double RandGauss::operator()() {
  return fire( defaultMean, defaultStdDev );
}

double RandGauss::operator()( double mean, double stdDev ) {
  return fire( mean, stdDev );
}

double RandGauss::shoot()
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( getFlag() ) {
    setFlag(false);
    double x = getVal();
    return x; 
    // return getVal();
  } 

  double r;
  double v1,v2,fac,val;
  HepRandomEngine* anEngine = HepRandom::getTheEngine();

  do {
    v1 = 2.0 * anEngine->flat() - 1.0;
    v2 = 2.0 * anEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt(-2.0*log(r)/r);
  val = v1*fac;
  setVal(val);
  setFlag(true);
  return v2*fac;
}

void RandGauss::shootArray( const int size, double* vect,
                            double mean, double stdDev )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(mean,stdDev);
}

double RandGauss::shoot( HepRandomEngine* anEngine )
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( getFlag() ) {
    setFlag(false);
    return getVal();
  }

  double r;
  double v1,v2,fac,val;

  do {
    v1 = 2.0 * anEngine->flat() - 1.0;
    v2 = 2.0 * anEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt( -2.0*log(r)/r);
  val = v1*fac;
  setVal(val);
  setFlag(true);
  return v2*fac;
}

void RandGauss::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double mean, double stdDev )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(anEngine,mean,stdDev);
}

double RandGauss::normal()
{
  // Gaussian random numbers are generated two at the time, so every other
  // time this is called we just return a number generated the time before.

  if ( set ) {
    set = false;
    return nextGauss;
  }

  double r;
  double v1,v2,fac,val;

  do {
    v1 = 2.0 * localEngine->flat() - 1.0;
    v2 = 2.0 * localEngine->flat() - 1.0;
    r = v1*v1 + v2*v2;
  } while ( r > 1.0 );

  fac = sqrt(-2.0*log(r)/r);
  val = v1*fac;
  nextGauss = val;
  set = true;
  return v2*fac;
}

void RandGauss::fireArray( const int size, double* vect)
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( defaultMean, defaultStdDev );
}

void RandGauss::fireArray( const int size, double* vect,
                           double mean, double stdDev )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire( mean, stdDev );
}

void RandGauss::saveEngineStatus ( const char filename[] ) {

  // First save the engine status just like the base class would do:
  getTheEngine()->saveStatus( filename );

  // Now append the cached variate, if any:

  std::ofstream outfile ( filename, std::ios::app );

  if ( getFlag() ) {
    std::vector<unsigned long> t(2);
    t = DoubConv::dto2longs(getVal());
    outfile << "RANDGAUSS CACHED_GAUSSIAN: Uvec " 
	    << getVal() << " " << t[0] << " " << t[1] << "\n";
  } else {
    outfile << "RANDGAUSS NO_CACHED_GAUSSIAN: 0 \n" ;
  }

} // saveEngineStatus

void RandGauss::restoreEngineStatus( const char filename[] ) {

  // First restore the engine status just like the base class would do:
  getTheEngine()->restoreStatus( filename );

  // Now find the line describing the cached variate:

  std::ifstream infile ( filename, std::ios::in );
  if (!infile) return;

  char inputword[] = "NO_KEYWORD    "; // leaves room for 14 characters plus \0
  while (true) {
    infile.width(13);
    infile >> inputword;
    if (strcmp(inputword,"RANDGAUSS")==0) break;
    if (infile.eof()) break;
	// If the file ends without the RANDGAUSS line, that means this
	// was a file produced by an earlier version of RandGauss.  We will
	// replicated the old behavior in that case:  set_st is cleared.
  }

  // Then read and use the caching info:

  if (strcmp(inputword,"RANDGAUSS")==0) {
    char setword[40];	// the longest, staticFirstUnusedBit: has length 21
    infile.width(39);
    infile >> setword;  // setword should be CACHED_GAUSSIAN:
    if (strcmp(setword,"CACHED_GAUSSIAN:") ==0) {
      if (possibleKeywordInput(infile, "Uvec", nextGauss_st)) {
        std::vector<unsigned long> t(2);
        infile >> nextGauss_st >> t[0] >> t[1]; 
        nextGauss_st = DoubConv::longs2double(t); 
      }
      // is >> nextGauss_st encompassed by possibleKeywordInput
      setFlag(true);
    } else {
      setFlag(false);
      infile >> nextGauss_st; // because a 0 will have been output
    }
  } else {
    setFlag(false);
  }

} // restoreEngineStatus

  // Save and restore to/from streams
  
std::ostream & RandGauss::put ( std::ostream & os ) const {
  os << name() << "\n";
  int prec = os.precision(20);
  std::vector<unsigned long> t(2);
  os << "Uvec\n";
  t = DoubConv::dto2longs(defaultMean);
  os << defaultMean << " " << t[0] << " " << t[1] << "\n";
  t = DoubConv::dto2longs(defaultStdDev);
  os << defaultStdDev << " " << t[0] << " " << t[1] << "\n";
  if ( set ) {
    t = DoubConv::dto2longs(nextGauss);
    os << "nextGauss " << nextGauss << " " << t[0] << " " << t[1] << "\n";
	#ifdef TRACE_IO
	std::cout << "put(): nextGauss = " << nextGauss << "\n";
	#endif
  } else {
	#ifdef TRACE_IO
	std::cout << "put(): No nextGauss \n";
	#endif
    os << "no_cached_nextGauss \n";
  }
  os.precision(prec);
  return os;
} // put   

std::istream & RandGauss::get ( std::istream & is ) {
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
  std::string c1;
  std::string c2;
  if (possibleKeywordInput(is, "Uvec", c1)) {
    std::vector<unsigned long> t(2);
    is >> defaultMean >> t[0] >> t[1]; defaultMean = DoubConv::longs2double(t); 
    is >> defaultStdDev>>t[0]>>t[1]; defaultStdDev = DoubConv::longs2double(t); 
    std::string ng;
    is >> ng;
    set = false;
	#ifdef TRACE_IO
	if (ng != "nextGauss") 
	std::cout << "get(): ng = " << ng << "\n";
	#endif	
    if (ng == "nextGauss") {
      is >> nextGauss >> t[0] >> t[1]; nextGauss = DoubConv::longs2double(t);
	#ifdef TRACE_IO
	std::cout << "get(): nextGauss read back as " << nextGauss << "\n";
	#endif	
      set = true;
    }
    return is;
  }
  // is >> c1 encompassed by possibleKeywordInput
  is >> defaultMean >> c2 >> defaultStdDev;
  if ( (!is) || (c1 != "Mean:") || (c2 != "Sigma:") ) {
    std::cerr << "i/o problem while expecting to read state of a "
    	      << name() << " distribution\n"
	      << "default mean and/or sigma could not be read\n";
    return is;
  }
  is >> c1 >> c2 >> nextGauss;
  if ( (!is) || (c1 != "RANDGAUSS") ) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Failure when reading caching state of RandGauss\n";
    return is;
  }
  if (c2 == "CACHED_GAUSSIAN:") {
    set = true;
  } else if (c2 == "NO_CACHED_GAUSSIAN:") {
    set = false;  
  } else {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Unexpected caching state keyword of RandGauss:" << c2
	      << "\nistream is left in the badbit state\n";
  } 
  return is;
} // get

  // Static save and restore to/from streams
  
std::ostream & RandGauss::saveDistState ( std::ostream & os ) {
  int prec = os.precision(20);
  std::vector<unsigned long> t(2);
  os << distributionName() << "\n";
  os << "Uvec\n";
  if ( getFlag() ) {
    t = DoubConv::dto2longs(getVal());
    os << "nextGauss_st " << getVal() << " " << t[0] << " " << t[1] << "\n";
  } else {
    os << "no_cached_nextGauss_st \n";
  }
  os.precision(prec);
  return os;
}    

std::istream & RandGauss::restoreDistState ( std::istream & is ) {
  std::string inName;
  is >> inName;
  if (inName != distributionName()) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Mismatch when expecting to read static state of a "
    	      << distributionName() << " distribution\n"
	      << "Name found was " << inName
	      << "\nistream is left in the badbit state\n";
    return is;
  }
  std::string c1;
  std::string c2;
  if (possibleKeywordInput(is, "Uvec", c1)) {
    std::vector<unsigned long> t(2);
    std::string ng;
    is >> ng;
    setFlag (false);
    if (ng == "nextGauss_st") {
      is >> nextGauss_st >> t[0] >> t[1]; 
      nextGauss_st = DoubConv::longs2double(t);
      setFlag (true);
    }
    return is;
  }
  // is >> c1 encompassed by possibleKeywordInput
  is >> c2 >> nextGauss_st;
  if ( (!is) || (c1 != "RANDGAUSS") ) {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Failure when reading caching state of static RandGauss\n";
    return is;
  }
  if (c2 == "CACHED_GAUSSIAN:") {
    setFlag(true);
  } else if (c2 == "NO_CACHED_GAUSSIAN:") {
    setFlag(false);  
  } else {
    is.clear(std::ios::badbit | is.rdstate());
    std::cerr << "Unexpected caching state keyword of static RandGauss:" << c2
	      << "\nistream is left in the badbit state\n";
  } 
  return is;
} 

std::ostream & RandGauss::saveFullState ( std::ostream & os ) {
  HepRandom::saveFullState(os);
  saveDistState(os);
  return os;
}
  
std::istream & RandGauss::restoreFullState ( std::istream & is ) {
  HepRandom::restoreFullState(is);
  restoreDistState(is);
  return is;
}

}  // namespace CLHEP

