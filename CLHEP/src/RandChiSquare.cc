// $Id: RandChiSquare.cc,v 1.1 2010/02/16 20:35:18 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- RandChiSquare ---
//                      class implementation file
// -----------------------------------------------------------------------

// =======================================================================
// John Marraffino - Created: 12th May 1998
// M Fischler     - put and get to/from streams 12/10/04
// M Fischler	      - put/get to/from streams uses pairs of ulongs when
//			+ storing doubles avoid problems with precision 
//			4/14/05
// =======================================================================

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/RandChiSquare.h"
#include "CLHEP/Random/DoubConv.hh"
#include <cmath>	// for log()

namespace CLHEP {

std::string RandChiSquare::name() const {return "RandChiSquare";}
HepRandomEngine & RandChiSquare::engine() {return *localEngine;}

RandChiSquare::~RandChiSquare() {
  if ( deleteEngine ) delete localEngine;
}

RandChiSquare::RandChiSquare(const RandChiSquare& right)
 : defaultA(right.defaultA)
{;}

double RandChiSquare::shoot( HepRandomEngine *anEngine,  double a ) {
  return genChiSquare( anEngine, a );
}

double RandChiSquare::shoot( double a ) {
  HepRandomEngine *anEngine = HepRandom::getTheEngine();
  return genChiSquare( anEngine, a );
}

double RandChiSquare::fire( double a ) {
  return genChiSquare( localEngine, a );
}

void RandChiSquare::shootArray( const int size, double* vect,
                            double a ) {
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(a);
}

void RandChiSquare::shootArray( HepRandomEngine* anEngine,
                            const int size, double* vect,
                            double a )
{
   int i;

   for (i=0; i<size; ++i)
     vect[i] = shoot(anEngine,a);
}

void RandChiSquare::fireArray( const int size, double* vect) {
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire(defaultA);
}

void RandChiSquare::fireArray( const int size, double* vect,
                           double a ) {
   int i;

   for (i=0; i<size; ++i)
     vect[i] = fire(a);
}

double RandChiSquare::genChiSquare( HepRandomEngine *anEngine,
                                       double a ) {
/******************************************************************
 *                                                                *
 *        Chi Distribution - Ratio of Uniforms  with shift        *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION :   - chru samples a random number from the Chi       *
 *                distribution with parameter  a > 1.             *
 * REFERENCE :  - J.F. Monahan (1987): An algorithm for           *
 *                generating chi random variables, ACM Trans.     *
 *                Math. Software 13, 168-172.                     *
 * SUBPROGRAM : - anEngine  ... pointer to a (0,1)-Uniform        *
 *                engine                                          *
 *                                                                *
 * Implemented by R. Kremer, 1990                                 *
 ******************************************************************/

 static double a_in = -1.0,b,vm,vp,vd;
 double u,v,z,zz,r;

// Check for invalid input value

 if( a < 1 )  return (-1.0);

 if (a == 1)
  {
   for(;;)
    {
     u = anEngine->flat();
     v = anEngine->flat() * 0.857763884960707;
     z = v / u;
     if (z < 0) continue;
     zz = z * z;
     r = 2.5 - zz;
     if (z < 0.0) r = r + zz * z / (3.0 * z);
     if (u < r * 0.3894003915) return(z*z);
     if (zz > (1.036961043 / u + 1.4)) continue;
     if (2 * log(u) < (- zz * 0.5 )) return(z*z);
     }
   }
 else
  {
   if (a != a_in)
    {
     b = sqrt(a - 1.0);
     vm = - 0.6065306597 * (1.0 - 0.25 / (b * b + 1.0));
     vm = (-b > vm)? -b : vm;
     vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
     vd = vp - vm;
     a_in = a;
     }
   for(;;)
    {
     u = anEngine->flat();
     v = anEngine->flat() * vd + vm;
     z = v / u;
     if (z < -b) continue;
     zz = z * z;
     r = 2.5 - zz;
     if (z < 0.0) r = r + zz * z / (3.0 * (z + b));
     if (u < r * 0.3894003915) return((z + b)*(z + b));
     if (zz > (1.036961043 / u + 1.4)) continue;
     if (2 * log(u) < (log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) return((z + b)*(z + b));
     }
   }
}

std::ostream & RandChiSquare::put ( std::ostream & os ) const {
  int pr=os.precision(20);
  std::vector<unsigned long> t(2);
  os << " " << name() << "\n";
  os << "Uvec" << "\n";
  t = DoubConv::dto2longs(defaultA);
  os << defaultA << " " << t[0] << " " << t[1] << "\n";
  os.precision(pr);
  return os;
#ifdef REMOVED
  int pr=os.precision(20);
  os << " " << name() << "\n";
  os << defaultA << "\n";
  os.precision(pr);
  return os;
#endif
}

std::istream & RandChiSquare::get ( std::istream & is ) {
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
    return is;
  }
  // is >> defaultA encompassed by possibleKeywordInput
  return is;
}



}  // namespace CLHEP

