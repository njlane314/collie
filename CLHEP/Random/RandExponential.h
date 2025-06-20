// $Id: RandExponential.h,v 1.1 2010/02/16 20:35:17 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- RandExponential ---
//                          class header file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// Class defining methods for shooting exponential distributed random
// values, given a mean (default mean = 1).
// Default mean is used for operator()().

// =======================================================================
// Gabriele Cosmo - Created: 5th September 1995
//                - Added methods to shoot arrays: 28th July 1997
// J.Marraffino   - Added default mean as attribute and
//                  operator() with mean: 16th Feb 1998
// M Fischler      - put and get to/from streams 12/10/04
// =======================================================================

#ifndef RandExponential_h
#define RandExponential_h 1

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/Random.h"

namespace CLHEP {

/**
 * @author <Gabriele.Cosmo@cern.ch>
 * @ingroup random
 */
class RandExponential : public HepRandom {

public:

  inline RandExponential ( HepRandomEngine& anEngine, double mean=1.0 );
  inline RandExponential ( HepRandomEngine* anEngine, double mean=1.0 );
  // These constructors should be used to instantiate a RandExponential
  // distribution object defining a local engine for it.
  // The static generator will be skipped using the non-static methods
  // defined below.
  // If the engine is passed by pointer the corresponding engine object
  // will be deleted by the RandExponential destructor.
  // If the engine is passed by reference the corresponding engine object
  // will not be deleted by the RandExponential destructor.

  virtual ~RandExponential();
  // Destructor

  // Static methods to shoot random values using the static generator

  static  double shoot();

  static  double shoot( double mean );

  static  void shootArray ( const int size, double* vect,
                            double mean=1.0 );

  //  Static methods to shoot random values using a given engine
  //  by-passing the static generator.

  static  inline double shoot( HepRandomEngine* anEngine );

  static  inline double shoot( HepRandomEngine* anEngine, double mean );

  static  void shootArray ( HepRandomEngine* anEngine, const int size,
                            double* vect, double mean=1.0 );

  //  Methods using the localEngine to shoot random values, by-passing
  //  the static generator.

  inline double fire();

  inline double fire( double mean );

  void fireArray ( const int size, double* vect );
  void fireArray ( const int size, double* vect, double mean );
  
  double operator()();
  double operator()( double mean );

  // Save and restore to/from streams
  
  std::ostream & put ( std::ostream & os ) const;
  std::istream & get ( std::istream & is );

  std::string name() const;
  HepRandomEngine & engine();

  static std::string distributionName() {return "RandExponential";}  
  // Provides the name of this distribution class
  
private:

  // Private copy constructor. Defining it here disallows use.
  RandExponential(const RandExponential& d);


  HepRandomEngine* localEngine;
  bool deleteEngine;
  double defaultMean;

};

}  // namespace CLHEP

#ifdef ENABLE_BACKWARDS_COMPATIBILITY
//  backwards compatibility will be enabled ONLY in CLHEP 1.9
using namespace CLHEP;
#endif

#include "CLHEP/Random/RandExponential.icc"

#endif
