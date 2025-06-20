// $Id: RanecuEngine.cc,v 1.1 2010/02/16 20:35:18 wfisher Exp $
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                        --- RanecuEngine ---
//                      class implementation file
// -----------------------------------------------------------------------
// This file is part of Geant4 (simulation toolkit for HEP).
//
// RANECU Random Engine - algorithm originally written in FORTRAN77
//                        as part of the MATHLIB HEP library.

// =======================================================================
// Gabriele Cosmo - Created - 2nd February 1996
//                - Minor corrections: 31st October 1996
//                - Added methods for engine status: 19th November 1996
//                - Added abs for setting seed index: 11th July 1997
//                - Modified setSeeds() to handle default index: 16th Oct 1997
//                - setSeed() now resets the engine status to the original
//                  values in the static table of HepRandom: 19th Mar 1998
// J.Marraffino   - Added stream operators and related constructor.
//                  Added automatic seed selection from seed table and
//                  engine counter: 16th Feb 1998
// Ken Smith      - Added conversion operators:  6th Aug 1998
// J. Marraffino  - Remove dependence on hepString class   13 May 1999
// M. Fischler    - Add endl to the end of saveStatus      10 Apr 2001
// M. Fischler    - In restore, checkFile for file not found    03 Dec 2004
// M. Fischler    - Methods for distrib. instance save/restore  12/8/04    
// M. Fischler    - split get() into tag validation and 
//                  getState() for anonymous restores           12/27/04    
// M. Fischler    - put/get for vectors of ulongs		3/14/05
// M. Fischler    - State-saving using only ints, for portability 4/12/05
//		    
// =======================================================================

#include "CLHEP/Random/defs.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Random/engineIDulong.h"
#include <string.h>
#include <cmath>
#include <stdlib.h>

namespace CLHEP {

static const int MarkerLen = 64; // Enough room to hold a begin or end marker. 

std::string RanecuEngine::name() const {return "RanecuEngine";}

// Number of instances with automatic seed selection
int RanecuEngine::numEngines = 0;

RanecuEngine::RanecuEngine()
: ecuyer_a(40014),ecuyer_b(53668),ecuyer_c(12211),
  ecuyer_d(40692),ecuyer_e(52774),ecuyer_f(3791),
  shift1(2147483563),shift2(2147483399),
  prec(4.6566128E-10 )
{
  int cycle = abs(int(numEngines/maxSeq));
  seq = abs(int(numEngines%maxSeq));
  numEngines += 1;
  theSeed = seq;
  long mask = ((cycle & 0x007fffff) << 8);
  for (int i=0; i<2; ++i) {
    for (int j=0; j<maxSeq; ++j) {
      HepRandom::getTheTableSeeds(table[j],j);
      table[j][i] ^= mask;
    }
  }
  theSeeds = &table[seq][0];
}

RanecuEngine::RanecuEngine(int index)
: ecuyer_a(40014),ecuyer_b(53668),ecuyer_c(12211),
  ecuyer_d(40692),ecuyer_e(52774),ecuyer_f(3791),
  shift1(2147483563),shift2(2147483399),
  prec(4.6566128E-10 )
{
  int cycle = abs(int(index/maxSeq));
  seq = abs(int(index%maxSeq));
  theSeed = seq;
  long mask = ((cycle & 0x000007ff) << 20);
  for (int j=0; j<maxSeq; ++j) {
    HepRandom::getTheTableSeeds(table[j],j);
    table[j][0] ^= mask;
    table[j][1] ^= mask;
  }
  theSeeds = &table[seq][0];
}

RanecuEngine::RanecuEngine(std::istream& is)
: ecuyer_a(40014),ecuyer_b(53668),ecuyer_c(12211),
  ecuyer_d(40692),ecuyer_e(52774),ecuyer_f(3791),
  shift1(2147483563),shift2(2147483399),
  prec(4.6566128E-10 )
{
   is >> *this;
}

RanecuEngine::~RanecuEngine() {}

RanecuEngine::RanecuEngine(const RanecuEngine &p)
: ecuyer_a(40014),ecuyer_b(53668),ecuyer_c(12211),
  ecuyer_d(40692),ecuyer_e(52774),ecuyer_f(3791),
  shift1(2147483563),shift2(2147483399),
  prec(4.6566128E-10 )
{
  if ((this != &p) && (&p)) {
    theSeed = p.getSeed();
    seq = p.seq;
    for (int i=0; i<2; ++i)
      for (int j=0; j<maxSeq; ++j)
        table[j][i] = p.table[j][i];
    theSeeds = &table[seq][0];
  }
}

RanecuEngine & RanecuEngine::operator = (const RanecuEngine &p)
{
  if ((this != &p) && (&p)) {
    theSeed = p.getSeed();
    seq = p.seq;
    for (int i=0; i<2; ++i)
      for (int j=0; j<maxSeq; ++j)
        table[j][i] = p.table[j][i];
    theSeeds = &table[seq][0];
  }
  return *this;
}

void RanecuEngine::setSeed(long index, int)
{
  seq = abs(int(index%maxSeq));
  theSeed = seq;
  HepRandom::getTheTableSeeds(table[seq],seq);
  theSeeds = &table[seq][0];
}

void RanecuEngine::setSeeds(const long* seeds, int pos)
{
  if (pos != -1) {
    seq = abs(int(pos%maxSeq));
    theSeed = seq;
  }
  // only positive seeds are allowed
  table[seq][0] = abs(seeds[0])%ecuyer_a;
  table[seq][1] = abs(seeds[1])%ecuyer_d;
  theSeeds = &table[seq][0];
}

void RanecuEngine::setIndex(long index)
{
  seq = abs(int(index%maxSeq));
  theSeed = seq;
  theSeeds = &table[seq][0];
}

void RanecuEngine::saveStatus( const char filename[] ) const
{
   std::ofstream outFile( filename, std::ios::out ) ;

  if (!outFile.bad()) {
    outFile << "Uvec\n";
    std::vector<unsigned long> v = put();
		     #ifdef TRACE_IO
			 std::cout << "Result of v = put() is:\n"; 
		     #endif
    for (unsigned int i=0; i<v.size(); ++i) {
      outFile << v[i] << "\n";
		     #ifdef TRACE_IO
			   std::cout << v[i] << " ";
			   if (i%6==0) std::cout << "\n";
		     #endif
    }
		     #ifdef TRACE_IO
			 std::cout << "\n";
		     #endif
  }
#ifdef REMOVED
   if (!outFile.bad()) {
     outFile << theSeed << std::endl;
     for (int i=0; i<2; ++i)
       outFile << table[theSeed][i] << " ";
     outFile << std::endl;
   }
#endif
}

void RanecuEngine::restoreStatus( const char filename[] )
{
  std::ifstream inFile( filename, std::ios::in);
  if (!checkFile ( inFile, filename, engineName(), "restoreStatus" )) {
    std::cerr << "  -- Engine state remains unchanged\n";
    return;
  }
  if ( possibleKeywordInput ( inFile, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long xin;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      inFile >> xin;
	       #ifdef TRACE_IO
	       std::cout << "ivec = " << ivec << "  xin = " << xin << "    ";
	       if (ivec%3 == 0) std::cout << "\n"; 
	       #endif
      if (!inFile) {
        inFile.clear(std::ios::badbit | inFile.rdstate());
        std::cerr << "\nJamesRandom state (vector) description improper."
	       << "\nrestoreStatus has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return;
      }
      v.push_back(xin);
    }
    getState(v);
    return;
  }

  if (!inFile.bad() && !inFile.eof()) {
//     inFile >> theSeed;  removed -- encompased by possibleKeywordInput
     for (int i=0; i<2; ++i)
       inFile >> table[theSeed][i];
     seq = int(theSeed);
  }
}

void RanecuEngine::showStatus() const
{
   std::cout << std::endl;
   std::cout << "--------- Ranecu engine status ---------" << std::endl;
   std::cout << " Initial seed (index) = " << theSeed << std::endl;
   std::cout << " Current couple of seeds = "
	     << table[theSeed][0] << ", "
	     << table[theSeed][1] << std::endl;
   std::cout << "----------------------------------------" << std::endl;
}

double RanecuEngine::flat()
{
   const int index = seq;
   long seed1 = table[index][0];
   long seed2 = table[index][1];

   int k1 = (int)(seed1/ecuyer_b);
   int k2 = (int)(seed2/ecuyer_e);

   seed1 = ecuyer_a*(seed1-k1*ecuyer_b)-k1*ecuyer_c;
   if (seed1 < 0) seed1 += shift1;
   seed2 = ecuyer_d*(seed2-k2*ecuyer_e)-k2*ecuyer_f;
   if (seed2 < 0) seed2 += shift2;

   table[index][0] = seed1;
   table[index][1] = seed2;

   long diff = seed1-seed2;

   if (diff <= 0) diff += (shift1-1);
   return (double)(diff*prec);
}

void RanecuEngine::flatArray(const int size, double* vect)
{
   const int index = seq;
   long seed1 = table[index][0];
   long seed2 = table[index][1];
   int k1, k2;
   register int i;

   for (i=0; i<size; ++i)
   {
     k1 = (int)(seed1/ecuyer_b);
     k2 = (int)(seed2/ecuyer_e);

     seed1 = ecuyer_a*(seed1-k1*ecuyer_b)-k1*ecuyer_c;
     if (seed1 < 0) seed1 += shift1;
     seed2 = ecuyer_d*(seed2-k2*ecuyer_e)-k2*ecuyer_f;
     if (seed2 < 0) seed2 += shift2;

     long diff = seed1-seed2;
     if (diff <= 0) diff += (shift1-1);

     vect[i] = (double)(diff*prec);
   }
   table[index][0] = seed1;
   table[index][1] = seed2;
}

RanecuEngine::operator unsigned int() {
   const int index = seq;
   long seed1 = table[index][0];
   long seed2 = table[index][1];

   int k1 = (int)(seed1/ecuyer_b);
   int k2 = (int)(seed2/ecuyer_e);

   seed1 = ecuyer_a*(seed1-k1*ecuyer_b)-k1*ecuyer_c;
   if (seed1 < 0) seed1 += shift1;
   seed2 = ecuyer_d*(seed2-k2*ecuyer_e)-k2*ecuyer_f;
   if (seed2 < 0) seed2 += shift2;

   table[index][0] = seed1;
   table[index][1] = seed2;
   long diff = seed1-seed2;
   if( diff <= 0 ) diff += (shift1-1);

   return ((diff << 1) | (seed1&1))& 0xffffffff;
}

std::ostream & RanecuEngine::put( std::ostream& os ) const
{
   char beginMarker[] = "RanecuEngine-begin";
  os << beginMarker << "\nUvec\n";
  std::vector<unsigned long> v = put();
  for (unsigned int i=0; i<v.size(); ++i) {
     os <<  v[i] <<  "\n";
  }
  return os;  
#ifdef REMOVED 
   char endMarker[]   = "RanecuEngine-end";
   os << " " << beginMarker << "\n";
   os << theSeed << " ";
   for (int i=0; i<2; ++i) {
     os << table[theSeed][i] << "\n";
   }
   os << endMarker << "\n";
   return os;
#endif
}

std::vector<unsigned long> RanecuEngine::put () const {
  std::vector<unsigned long> v;
  v.push_back (engineIDulong<RanecuEngine>());
  v.push_back(static_cast<unsigned long>(theSeed));
  v.push_back(static_cast<unsigned long>(table[theSeed][0]));
  v.push_back(static_cast<unsigned long>(table[theSeed][1]));
  return v;
}

std::istream & RanecuEngine::get ( std::istream& is )
{
  char beginMarker [MarkerLen];

  is >> std::ws;
  is.width(MarkerLen);  // causes the next read to the char* to be <=
			// that many bytes, INCLUDING A TERMINATION \0 
			// (Stroustrup, section 21.3.2)
  is >> beginMarker;
  if (strcmp(beginMarker,"RanecuEngine-begin")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nInput stream mispositioned or"
	       << "\nRanecuEngine state description missing or"
	       << "\nwrong engine type found." << std::endl;
     return is;
   }
  return getState(is);
}

std::string RanecuEngine::beginTag ( )  { 
  return "RanecuEngine-begin"; 
}

std::istream & RanecuEngine::getState ( std::istream& is )
{
  if ( possibleKeywordInput ( is, "Uvec", theSeed ) ) {
    std::vector<unsigned long> v;
    unsigned long uu;
    for (unsigned int ivec=0; ivec < VECTOR_STATE_SIZE; ++ivec) {
      is >> uu;
      if (!is) {
        is.clear(std::ios::badbit | is.rdstate());
        std::cerr << "\nRanecuEngine state (vector) description improper."
		<< "\ngetState() has failed."
	       << "\nInput stream is probably mispositioned now." << std::endl;
        return is;
      }
      v.push_back(uu);
    }
    getState(v);
    return (is);
  }

//  is >> theSeed;  Removed, encompassed by possibleKeywordInput()
  char endMarker   [MarkerLen];
   for (int i=0; i<2; ++i) {
     is >> table[theSeed][i];
   }
  is >> std::ws;
  is.width(MarkerLen);  
  is >> endMarker;
  if (strcmp(endMarker,"RanecuEngine-end")) {
     is.clear(std::ios::badbit | is.rdstate());
     std::cerr << "\nRanecuEngine state description incomplete."
	       << "\nInput stream is probably mispositioned now." << std::endl;
     return is;
   }

   seq = int(theSeed);
   return is;
}

bool RanecuEngine::get (const std::vector<unsigned long> & v) {
  if ((v[0] & 0xffffffffUL) != engineIDulong<RanecuEngine>()) {
    std::cerr << 
    	"\nRanecuEngine get:state vector has wrong ID word - state unchanged\n";
    return false;
  }
  return getState(v);
}

bool RanecuEngine::getState (const std::vector<unsigned long> & v) {
  if (v.size() != VECTOR_STATE_SIZE ) {
    std::cerr << 
    	"\nRanecuEngine get:state vector has wrong length - state unchanged\n";
    return false;
  }
  theSeed           = v[1];
  table[theSeed][0] = v[2];
  table[theSeed][1] = v[3];
  seq = int(theSeed);
  return true;
}


}  // namespace CLHEP
