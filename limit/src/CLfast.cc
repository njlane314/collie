#include <math.h>
#include <algorithm>
#include <sys/time.h>
#include "CLfast.hh"
#include "timeBasedSeed.hh" // m. fischler 1/8/09

template <class T>
inline const T& MAXof(const T& a, const T& b) { return (a>b)?a:b; }

CLfast::CLfast() {
  timeval a;
  gettimeofday(&a,NULL); 
  rp_ = new RandPoisson(new MTwistEngine(timeBasedSeed()),1);
  // was: rp_ = new RandPoisson(new MTwistEngine((int)(a.tv_usec%100000+12)),1);
  // in 1.20:  rp_ = new RandPoisson(new MTwistEngine(a.tv_usec%3+2),1);
  iterations_ = 100000;
  nmcdone_ = 0;
  useHistoStats_ = false;
  medianExpected_ = true;
  setNoviceFlag(true);
}

bool CLfast::calculateCLs(const SigBkgdDist& sbd, CLpoint& CLs, int effort) {

  if (effort==LEVEL_TEST) return doCLs(sbd,CLs,100);
  else if (effort==LEVEL_VERYVERYFAST) return doCLs(sbd,CLs,5000);
  else if (effort==LEVEL_VERYFAST) return doCLs(sbd,CLs,15000);
  else if (effort==LEVEL_FAST) return doCLs(sbd,CLs,25000);
  else if (effort==LEVEL_STANDARD) return doCLs(sbd,CLs,50000);
  else if (effort==LEVEL_FINE) return doCLs(sbd,CLs,100000);
  else if (effort==LEVEL_VERYFINE) return doCLs(sbd,CLs,200000);
  else if (effort==LEVEL_VERYVERYFINE) return doCLs(sbd,CLs,500000);
  else // LEVEL_FAST
    return doCLs(sbd,CLs,25000);
}


bool CLfast::doCLs(const SigBkgdDist& sbd, CLpoint& CLs, int its) {

  double llr_d=0, llr_b_expected=0, llr_sb_expected=0, llr_sb=0, llr_b=0;
  SigBkgdDist sblike(sbd), blike(sbd);

  const double* xback=sblike.bkgd();
  const double* xsig=sblike.signal();
  
  int nbins=sblike.nbins();
  int numer=0; int denom=0;
  int denom_bexp=0;
  
  //cache the data LLR value before resetting the data arrays
  llr_d=sblike.calculateLLR();
  
  //adjust the data predictions for the two hypotheses
  for(int i=0; i<nbins; i++) {
    sblike.data(i)= xback[i]+xsig[i];
    blike.data(i)= xback[i];
  }

  llr_sb_expected=sblike.calculateLLR();
  llr_b_expected=blike.calculateLLR();

  CLs.llrobs=llr_d;
  CLs.llrb  =llr_b_expected;
  CLs.llrsb =llr_sb_expected;


  int nmc;
  double zback, zsig;
  for (nmc=0; nmc<2*its && (denom<its/2.0 || nmc<its); nmc++) {
    for (int j=0; j<nbins; j++) {
      sblike.data(j)=0; blike.data(j)=0;
      if(xback[j]==0 && xsig[j]==0) continue;   //skip zero bins...
      if (xback[j]>0) {
	zback=rp_->fire(xback[j]);
      } else zback=0;
      if (xsig[j]>0) {
	zsig = rp_->fire(xsig[j]);
      } else zsig=0;

      blike.data(j)=zback;
      sblike.data(j)=zsig+zback;
    }

    llr_b = blike.calculateLLR();
    llr_sb=sblike.calculateLLR();

    if (llr_d<=llr_b){   denom++;  // background "less than" data
      if (llr_d<=llr_sb) numer++; // s+b "less than" data
    }
    if (llr_b_expected<=llr_b) denom_bexp++;// background "less than" med expected b
    if (nmc<MaxKeep_) {
       btrial_[nmc]=llr_b;
      sbtrial_[nmc]=llr_sb;
    }
  }


  nmcdone_=nmc;
  CLs.clb_obs  =(double(denom)/double(nmcdone_));
  CLs.clsb_obs =(double(numer)/double(nmcdone_));
  CLs.clb_med  =(double(denom_bexp)/double(nmcdone_));

  double cls_obs=1.0;
  if (denom>0) cls_obs=1.0-double(numer)/double(denom); // fraction of trials where data < s+b
  if (sbd.totSignal()<1e-5) cls_obs=0.0;
  CLs.cls_obs=cls_obs;
  
  int lim=(nmcdone_<MaxKeep_)?nmcdone_:MaxKeep_;
  
  CLs.nTrials_obs = nmcdone_;
  CLs.nTrials_exp = lim;

  std::sort(btrial_,btrial_+lim);
  std::sort(sbtrial_,sbtrial_+lim);
  
  CLs.cls_med=-1.0;
  CLs.cls_med_m1s=-1.0;
  CLs.cls_med_m2s=-1.0;
  CLs.cls_med_p1s=-1.0;
  CLs.cls_med_p2s=-1.0;
  
  int m1s_pos = int(lim*(1.0-0.3173/2.0)+0.5)+1;
  int p1s_pos = int(lim*(    0.3173/2.0)+0.5)+1;
  int m2s_pos = int(lim*(1.0-0.0455/2.0)+0.5)+1;
  int p2s_pos = int(lim*(    0.0455/2.0)+0.5)+1;

  CLs.llrsb_m2s=sbtrial_[m2s_pos];
  CLs.llrsb_m1s=sbtrial_[m1s_pos];
  CLs.llrsb_p1s=sbtrial_[p1s_pos];
  CLs.llrsb_p2s=sbtrial_[p2s_pos];

  CLs.llrb_m2s=btrial_[m2s_pos];
  CLs.llrb_m1s=btrial_[m1s_pos];
  CLs.llrb_p1s=btrial_[p1s_pos];
  CLs.llrb_p2s=btrial_[p2s_pos];
  
  double val=btrial_[lim/2];
  for (int i=0; i<lim; i++) {
    if (sbtrial_[i]>=val && CLs.cls_med<0) { 
      CLs.clsb_med = double(lim-i)/double(lim); 
       CLs.cls_med=MAXof(0.0,double(i)/(lim/2)-1.0); 
    }
    if (sbtrial_[i]>=btrial_[m1s_pos] && CLs.cls_med_m1s<0){ 
       CLs.cls_med_m1s=1.0-double(lim-i)/double(lim-m1s_pos);
      CLs.clsb_med_m1s=double(lim-i)/double(lim); 
    }
    if (sbtrial_[i]>=btrial_[m2s_pos] && CLs.cls_med_m2s<0){ 
       CLs.cls_med_m2s=1.0-double(lim-i)/double(lim-m2s_pos);
      CLs.clsb_med_m2s=double(lim-i)/double(lim);
    }
    if (sbtrial_[i]>=btrial_[p1s_pos] && CLs.cls_med_p1s<0){ 
       CLs.cls_med_p1s=1.0-double(lim-i)/double(lim-p1s_pos);
      CLs.clsb_med_p1s=double(lim-i)/double(lim);
    }
    if (sbtrial_[i]>=btrial_[p2s_pos] && CLs.cls_med_p2s<0){ 
       CLs.cls_med_p2s=1.0-double(lim-i)/double(lim-p2s_pos);
      CLs.clsb_med_p2s=double(lim-i)/double(lim);
    }
  }
  if (CLs.cls_med==-1.0) CLs.cls_med=1.0;
  if (CLs.cls_med_m1s==-1.0) CLs.cls_med_m1s=1.0;
  if (CLs.cls_med_m2s==-1.0) CLs.cls_med_m2s=1.0;
  if (CLs.cls_med_p1s==-1.0) CLs.cls_med_p1s=1.0;
  if (CLs.cls_med_p2s==-1.0) CLs.cls_med_p2s=1.0;
  
  if (CLs.cls_med>1.0 && CLs.cls_med<1.001) CLs.cls_med=1.0;
  if (CLs.cls_med<0.0) CLs.cls_med=0.0;  
  CLs.clsb_med = (1.0-CLs.cls_med)/2.0;

  
  if(!medianExpected_){
    val = llr_b_expected;
    CLs.cls_med = -1.0;
    for (int i=0; i<lim; i++) {
      if (sbtrial_[i]>=val && CLs.cls_med<0){
	CLs.clsb_med = 1.0-double(i)/double(lim);   
	CLs.cls_med=MAXof(0.0,1-CLs.clsb_med/CLs.clb_med);	
      }
    }
  }

  CLs.clb_sb=-1.0;
  CLs.clb_sb_p1s=-1.0;
  CLs.clb_sb_p2s=-1.0;
  CLs.clb_sb_m1s=-1.0;
  CLs.clb_sb_m2s=-1.0;
  for (int i=0; i<lim; i++) {    
    if (btrial_[i]>CLs.llrsb && CLs.clb_sb<0)         CLs.clb_sb    =1.0-double(i)/double(lim); 
    if (btrial_[i]>CLs.llrsb_p1s && CLs.clb_sb_p1s<0) CLs.clb_sb_p1s=1.0-double(i)/double(lim); 
    if (btrial_[i]>CLs.llrsb_p2s && CLs.clb_sb_p2s<0) CLs.clb_sb_p2s=1.0-double(i)/double(lim); 
    if (btrial_[i]>CLs.llrsb_m1s && CLs.clb_sb_m1s<0) CLs.clb_sb_m1s=1.0-double(i)/double(lim); 
    if (btrial_[i]>CLs.llrsb_m2s && CLs.clb_sb_m2s<0) CLs.clb_sb_m2s=1.0-double(i)/double(lim); 
  }
  if (CLs.clb_sb==-1.0) CLs.clb_sb=0;
  if (CLs.clb_sb_m1s==-1.0) CLs.clb_sb_m1s=0;
  if (CLs.clb_sb_m2s==-1.0) CLs.clb_sb_m2s=0;
  if (CLs.clb_sb_p1s==-1.0) CLs.clb_sb_p1s=0;
  if (CLs.clb_sb_p2s==-1.0) CLs.clb_sb_p2s=0;

  return true;
}

void CLfast::configure(const char* options) {
  if (strstr(options,"its=")!=NULL) {
    iterations_=atoi(strstr(options,"its=")+4);
  }
}

TH1D* CLfast::getLLRdist_sb(const char* title,int bins,double min,double max){
  
  if(max < min) {
    min = sbtrial_[0];
    max = sbtrial_[(nmcdone_<MaxKeep_)?nmcdone_-1:MaxKeep_-1];
    std::cout << " SB MIN: " << min << " MAX " << max << std::endl;
  }
  TH1D* out = new TH1D(title,title,bins,min,max);
  int lim=(nmcdone_<MaxKeep_)?nmcdone_:MaxKeep_;

  for(int i=0; i<lim; i++)
    out->Fill(sbtrial_[i]);
  
  return out;
}

TH1D* CLfast::getLLRdist_b(const char* title,int bins,double min,double max){
  
  if(max < min) {
    min = btrial_[0];
    max = btrial_[(nmcdone_<MaxKeep_)?nmcdone_-1:MaxKeep_-1];
    std::cout << " B MIN: " << min << " MAX " << max << std::endl;
  }
  TH1D* out = new TH1D(title,title,bins,min,max);
  int lim=(nmcdone_<MaxKeep_)?nmcdone_:MaxKeep_;

  for(int i=0; i<lim; i++)
    out->Fill(btrial_[i]);
  
  return out;
}

void CLfast::noviceWarning(){  
  if(useHistoStats_){
    printf("\n*****************************************************************\n");
    printf("CLfast Warning:  You're operating in novice mode and you're using\nstatistical uncertainties from your input histograms.  Are you\nsure you've set the values properly in the histograms?\n");
    printf("*****************************************************************\n\n");
  }
  else{
    printf("\n*****************************************************************\n");
    printf("CLfast Warning:  The histogram statistical uncertainty flag is turned\noff by default.  This calculation does not include statistical uncertainties.\n");
    printf("*****************************************************************\n\n");  
  }
}
