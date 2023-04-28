#ifndef CLfast_hh_included
#define CLfast_hh_included

#include "SigBkgdDist.hh"
#include "CLcompute.hh"

/** Confidence level calculation by Monte Carlo trials. */
/**
   This class calculates confidence levels via integration of Poisson
   log-likelihood ratios for background-only (NULL) hypotheses and
   signal plus background (TEST) hypotheses.  In this class,
   uncertainties on nuisance parameters are **NOT** taken into
   account.  This class is intended to be a fast check on the validity
   of the inputs and the relative signal/background/data agreement.
   Do not use this class to generate limits for public consumption.
**/

class CLfast : public CLcompute {
public:
  /// constructor
  CLfast();
  
  /** Defined options include:
      "its=(number of Monte Carlo trials to throw)"
  */
  virtual void configure(const char* options);
  virtual bool calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort);
  
  // Calculate statistical uncertainty using contents of histograms.  Default = false
  //    Only turn this on if you're sure your histograms contain the correct
  //    statistical errors per bin!!
  void useHistoStats(bool use) { useHistoStats_ = use; }

  // Get the total number of MC iterations performed
  inline int getNMCiter() const { return nmcdone_; }

  TH1D* getLLRdist_sb(const char* title,int bins,double min,double max);
  TH1D* getLLRdist_b(const char* title,int bins,double min,double max);

  void noviceWarning();

private:
  bool doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its);
  
  static const int MaxKeep_ = (4*1024*1024)/sizeof(double)/2; // allocate 5 MB to this task 
  
  double btrial_[MaxKeep_];
  double sbtrial_[MaxKeep_];
  
  RandPoisson* rp_;
  
  int iterations_;
  
  int nmcdone_;
  
  bool useHistoStats_;  
};

#endif
