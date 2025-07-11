// $Id: RandGauss.h,v 1.1 2010/02/16 20:35:17 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandGauss ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).

// Class defining methods for shooting gaussian distributed random values,
// given a mean (default=0) or specifying also a deviation (default=1).
// Gaussian random numbers are generated two at the time, so every
// other time shoot is called the number returned is the one generated the
// time before.
// Default values are used for operator()().

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Minor corrections: 31st October 1996
//                - Added methods to shoot arrays: 28th July 1997
// J.Marraffino   - Added default arguments as attributes and
//                  operator() with arguments. Introduced method normal()
//                  for computation in fire(): 16th Feb 1998
// Gabriele Cosmo - Relocated static data from HepRandom: 5th Jan 1999
// M Fischler     - put and get to/from streams 12/8/04
// =======================================================================

#ifndef RandGauss_h
#define RandGauss_h 1

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/Random.h"

namespace CLHEP {

/**
 * @author
 * @ingroup random
 */
class RandGauss : public HepRandom {

public:

  inline RandGauss ( HepRandomEngine& anEngine, double mean=0.0,
                                                double stdDev=1.0 );
  inline RandGauss ( HepRandomEngine* anEngine, double mean=0.0,
                                                double stdDev=1.0 );
  // These constructors should be used to instantiate a RandGauss
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandGauss destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandGauss destructor.

  virtual ~RandGauss();
  // Destructor

  // Static methods to shoot random values using the static generator

  static  double shoot();

  static  inline double shoot( double mean, double stdDev );

  static  void shootArray ( const int size, double* vect,
                            double mean=0.0, double stdDev=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  double shoot( HepRandomEngine* anEngine );

  static  inline double shoot( HepRandomEngine* anEngine, 
                                  double mean, double stdDev );

  static  void shootArray ( HepRandomEngine* anEngine, const int size,
                            double* vect, double mean=0.0,
                            double stdDev=1.0 );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  double fire();

  inline double fire( double mean, double stdDev );
  
  void fireArray ( const int size, double* vect);
  void fireArray ( const int size, double* vect,
                   double mean, double stdDev );

  virtual double operator()();
  virtual double operator()( double mean, double stdDev );

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandGauss";}  
  // Provides the name of this distribution class
    
  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );
  
  //  Methods setFlag(false) and setF(false) if invoked in the client
  //  code before shoot/fire will force generation of a new couple of
  //  values.

  static  bool getFlag() {return set_st;}

  static  void setFlag( bool val ) {set_st = val;}

  bool getF() const {return set;}
  
  void setF( bool val ) {set = val;}

  // Methods overriding the base class static saveEngineStatus ones,
  // by adding extra data so that save in one program, then further gaussians,
  // will produce the identical sequence to restore in another program, then 
  // generating gaussian randoms there 

  static void saveEngineStatus( const char filename[] = "Config.conf" );
  // Saves to file the current status of the current engine.

  static void restoreEngineStatus( const char filename[] = "Config.conf" );
  // Restores a saved status (if any) for the current engine.

  static std::ostream& saveFullState ( std::ostream & os );
  // Saves to stream the state of the engine and cached data.

  static std::istream& restoreFullState ( std::istream & is );
  // Restores from stream the state of the engine and cached data.

  static std::ostream& saveDistState ( std::ostream & os );
  // Saves to stream the state of the cached data.

  static std::istream& restoreDistState ( std::istream & is );
  // Restores from stream the state of the cached data.


protected:

  // Protected copy constructor. Defining it here disallows user use.
  RandGauss(const RandGauss& d);

  static  double getVal() {return nextGauss_st;}

  static  void setVal( double nextVal ) {nextGauss_st = nextVal;}

  double normal();

  double defaultMean;
  double defaultStdDev;

  HepRandomEngine* localEngine;

private:

  bool deleteEngine, set;
  double nextGauss;

  // static data
  static bool set_st;
  static double nextGauss_st;

};

}  // namespace CLHEP

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
using namespace CLHEP;
#endif

#include "CLHEP/Random/RandGauss.icc"

#endif
