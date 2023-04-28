#include "ProfileLH.hh"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "cfortran/cfortran.h"
#include "cfortran/minuit.h"
#include "TrialPoint.hh"

#ifndef ProfileParams_DEF
#define ProfileParams_DEF

SigBkgdDist* pf_sbd;
std::vector<double> pf_params; 
std::vector<uint> pf_index;
uint pf_nParm;
uint pf_iter;
uint pf_idx;
bool pf_fitSig;
bool pf_sigFloat;
double pf_sigLLR;
double pf_scl;
TrialPoint pf_lastMinuitTrialPoint;
#endif // ProfileParams_DEF

///Standard constructor
ProfileLH::ProfileLH(){
  setup();
}

///Constructor with SBD initializer
ProfileLH::ProfileLH(SigBkgdDist* asbd){
  if(asbd==NULL){
    printf("ProfileLH::ProfileLH, the input SBD is NULL!\n");
    return;
  }
  setup();
  setModel(asbd);
}



void ProfileLH::setup(){
  pf_sbd       = NULL;
  pf_params.clear();  // USED TO BE:  pf_params = NULL;
  pf_index.clear();
  pf_paramErrs.clear();
  pf_nParm = 0;
  pf_iter = 0;
  pf_fitSig   = true;
  pf_sigFloat = false;  
  pf_sigLLR = 1e6;
  pf_scl = 1.0;

  m_init      = false;
  m_fitTest   = false;
  m_fitSig    = true;
  m_verbose   = false;
  m_sigFloat  = false;  
  m_sigScale    = 1.0;  
  m_sigScaleNom = 1.0;
  m_sigScaleErr = 0.0;
  m_sigScaleErrP   = 0.0;
  m_sigScaleErrM   = 0.0;
  m_sigScaleParab  = 0.0;
  m_sigScaleGlobCC = 0.0;
  m_sigLLR      = 1e6;
  m_chi2min     = -1.0;
  m_edm         = -1.0;
  m_lun         = 99;
  m_status      = -1;
  m_maxIter     = 250000;

  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNames.clear();
  m_systValues.clear();
  m_systErrors.clear();
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  m_condI.sig=0;
  m_condI.bkgd =0;
  m_condI.data =0;
  m_condI.chi2 =0;
  m_condI.chi2_bins =0;
  m_condI.chi2_syst =0;

  m_condF.sig=0;
  m_condF.bkgd =0;
  m_condF.data =0;
  m_condF.chi2 =0;
  m_condF.chi2_bins =0;
  m_condF.chi2_syst =0;
}

void ProfileLH::setModel(SigBkgdDist* asbd){  
  pf_sbd = asbd;
  pf_fitSig = m_fitSig;
  pf_sigLLR = m_sigLLR;
  fillSyst();
}

int ProfileLH::getNiterations(){
  return pf_iter;
}

void ProfileLH::getMNerrs(){
  if(!m_init) { printf("ProfileLH::getMNerrs, not initialized\n"); return;}
  
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();

  //Deriving MINOS errors... 
  double arglist[2];
  arglist[0] = m_maxIter;
  arglist[1] = 1.0;
  int narg   = 2;
  int ierflg = 0;
  char* title = "MINOS";
  
  MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);

  double errP, errM, eparab, globcc;
  for(uint num=1;num<=pf_nParm; ++num){
    arglist[1] = num;
    MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);
    MNERRS(num,errP,errM,eparab,globcc);
    m_systErrorsUp.push_back(errP);
    m_systErrorsDn.push_back(errM);
    m_systErrorsParab.push_back(eparab);
    m_systErrorsGlobCC.push_back(globcc);
  }

  if(m_sigFloat){
    arglist[1] = pf_nParm+1;
    MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);
    MNERRS(pf_nParm+1,errP,errM,eparab,globcc);
    //    m_sigScaleErr = (fabs(errP)+fabs(errM))/2.0;
    m_sigScaleErrP   = errP;
    m_sigScaleErrM   = errM;
    m_sigScaleParab  = eparab;
    m_sigScaleGlobCC = globcc;
  }

  return;
}

bool ProfileLH::getEMat(double* emat){
  
  if(!m_init) { printf("ProfileLH::getEMat, not initialized\n"); return false;}
  
  m_systErrorsUp.clear();
  m_systErrorsDn.clear();
  m_systErrorsParab.clear();
  m_systErrorsGlobCC.clear();
  
  //Deriving MINOS errors... 
  double arglist[2];
  arglist[0] = m_maxIter;
  arglist[1] = 1.0;
  char* title = "MINOS";
  int narg = 2;
  int ierflg=0;  
  MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);
  double par3;
  int par4, par5;
  int status;
  double chi2min;
  MNSTAT(chi2min,m_edm,par3,par4,par5,status);
  if(status!=3) printf("\nProfileLH::getEMat, Status = %d\n\n",status);

  MNEMAT(*emat,pf_nParm);  

  double errP, errM, eparab, globcc;
  for(uint num=1;num<=pf_nParm; ++num){
    arglist[1] = num;
    MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);
    MNERRS(num,errP,errM,eparab,globcc);
    m_systErrorsUp.push_back(errP);
    m_systErrorsDn.push_back(errM);
    m_systErrorsParab.push_back(eparab);
    m_systErrorsGlobCC.push_back(globcc);
  }

  if(m_sigFloat){
    arglist[1] = pf_nParm+1;
    MNEXCM(chiFun,title,arglist,narg,ierflg,NULL);
    MNERRS(pf_nParm+1,errP,errM,eparab,globcc);
    //    printf("mnerrs: %f, %f, %f, %f, %f\n",m_sigScaleErr,errP,errM,eparab,globcc);
    //    m_sigScaleErr = (fabs(errP)+fabs(errM))/2.0;    
    m_sigScaleErrP   = errP;
    m_sigScaleErrM   = errM;
    m_sigScaleParab  = eparab;
    m_sigScaleGlobCC = globcc;
  }

  return true;
}

void ProfileLH::fillSyst(){
  m_systNames.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();

  if(pf_sbd==NULL){
    printf("ProfileLH::fillSyst, You're trying to fit a profile with a NULL signal/background distribution!!\n"); 
    return;
  }
  
  //Make sure all the systematics are properly set
  if(!pf_sbd->isInit()){
    pf_sbd->varySystematics();    
    pf_sbd->setBaselineModel();
  }
  
  vector<string> chanNames = pf_sbd->getChannelNames();
  const CollieDistribution* dist;
  string systName ="";
  map<string,int> m_systMap;

  for(vector<string>::iterator its = chanNames.begin(); its!=chanNames.end(); ++its){
    for(uint b = 0; b<pf_sbd->getNbkgdDist(*its); ++b){
      
      dist = pf_sbd->getBkgdDist(*its,b);
      if(dist!=NULL){
	
	for(int s = 0; s<dist->getNsystematics(); ++s){
	  systName = dist->getSystName(s);
	  if(m_systMap.find(systName)!=m_systMap.end()) continue;
	  else{
	    m_systNames.push_back(systName);	    
	    if(pf_sbd->getSystFitFlag(systName)){
	      m_systNamesFitMap[systName] = pf_sbd->getSystIndex(systName);
	      m_systNamesFitBkgdMap[systName] = pf_sbd->getSystIndex(systName);
	      m_systNamesFit.push_back(systName);
	      m_systNamesFitBkgd.push_back(systName);
	    }
	    else{
	      if(m_verbose) printf("ProfileLH:: Not Fitting %s\n",dist->getSystName(s).c_str());
	    }
	    m_systMap[systName] = pf_sbd->getSystIndex(systName);
	  }
	}
      }
    }
  }

  for(vector<string>::iterator its = chanNames.begin(); its!=chanNames.end(); ++its){
    for(uint b = 0; b<pf_sbd->getNsigDist(*its); ++b){      
      //      dist = pf_sbd->getSigDist(*its,pf_sbd->getSigName(b));
      dist = pf_sbd->getSigDist(*its,b);
      if(dist!=NULL){	
	for(int s = 0; s<dist->getNsystematics(); ++s){	    
	  systName = dist->getSystName(s);
	  
	  if(m_systMap.find(systName)!=m_systMap.end()) continue;
	  else{
	    m_systNames.push_back(systName);
	    if(pf_sbd->getSystFitFlag(systName) && pf_fitSig){
	      m_systNamesFitMap[systName] = pf_sbd->getSystIndex(systName);
	      m_systNamesFit.push_back(systName);
	    }
	    else{
	      if(m_verbose) printf("ProfileLH:: Not Fitting %s\n",systName.c_str());
	    }
	    m_systMap[systName] = pf_sbd->getSystIndex(systName);
	  }
	}
      }// else printf("This sig dist is NULL???\n");
    }
  }

  //our adjustable copy of the syst params
  pf_params = std::vector<double>(m_systNames.size(),0.0); 
  pf_paramErrs = std::vector<double>(m_systNames.size(),0.0); 
  //  m_systErrors.clear();
  //  m_systValues.clear();  

  //collect existing parameters
  for(unsigned int i=0; i<m_systNames.size(); ++i)  pf_params[i] = pf_sbd->getSystFluctValue(i);
  
  initFitParams();

  return;
}

void ProfileLH::initFitParams(){
  //Make a storage container with only the values we're fitting...
  //   either background-only or sig+bkgd
  //   for bkgd-only, make index so we know where it fits in master copy.
  //   Allow individual fit parameters to be excluded by the user.
  
  pf_index.clear();
  map<string,int>::iterator im;
  
  if(m_fitSig){
    pf_nParm = m_systNamesFit.size();
    pf_index = std::vector<uint>(pf_nParm,0); 
    
    //now fill in the syst indices
    for(uint s=0; s<m_systNamesFit.size(); ++s){
      im = m_systNamesFitMap.find(m_systNamesFit[s]);
      pf_index[s] = im->second;
    }
  }
  else{
    pf_nParm = m_systNamesFitBkgd.size();
    pf_index = std::vector<uint>(pf_nParm,0); 
    
    //now fill in the syst indices
    for(uint s=0; s<m_systNamesFitBkgd.size(); ++s){
      im = m_systNamesFitBkgdMap.find(m_systNamesFitBkgd[s]);
      pf_index[s] = im->second;
    }
  }

  return;
}


bool ProfileLH::getMinosError(uint idx, double& errP, double& errM, double& parab, double& globCC){
  if(m_systErrorsUp.size()==0) getMNerrs();
  if(m_systErrorsUp.size()==0){ printf("ProfileLH::getMinosError, getMinosError failed...\n"); return false;}  //failed..
  if(idx > pf_nParm) return false; //failed..
  
  errP = m_systErrorsUp[idx];
  errM = m_systErrorsDn[idx];
  parab = m_systErrorsParab[idx];
  globCC = m_systErrorsGlobCC[idx];
  return true;
}

bool ProfileLH::getMinosErrXsec(double& errP, double& errM, double& parab, double& globCC){
  if(m_systErrorsUp.size()==0) getMNerrs();
  if(m_systErrorsUp.size()==0){ printf("ProfileLH::getMinosErrXsec, getMinosErrXsec failed...\n"); return false;}  //failed..

  errP = m_sigScaleErrP;
  errM = m_sigScaleErrM;
  parab = m_sigScaleParab;
  globCC = m_sigScaleGlobCC;
  return true;
}


double ProfileLH::getParamVal(uint p){
  if(p<0 || p>m_systNames.size()) return 0; 
  return pf_params[p];
}

double ProfileLH::getParamError(uint p){
  if(p<0 || p>pf_paramErrs.size()) return 0; 
  return pf_paramErrs[p];
}

int ProfileLH::getFitSystIndex(string name){
  map<string,int>::iterator iter;
  if(m_fitSig){
    for(iter = m_systNamesFitMap.begin(); iter!=m_systNamesFitMap.end(); ++iter)
      if(iter->first == name) return iter->second;
  }
  else{
    for(iter = m_systNamesFitBkgdMap.begin(); iter!=m_systNamesFitBkgdMap.end(); ++iter)
      if(iter->first == name) return iter->second;
  }
  return -1;
}

double ProfileLH::getFitParamErr(string name) const {
  
  if(m_fitSig){
    for(uint s=0; s<m_systNamesFit.size(); ++s){
      if(m_systNamesFit[s]==name) return m_systErrors[s];
    }
  }
  else{
    for(uint s=0; s<m_systNamesFitBkgd.size(); ++s){
      if(m_systNamesFitBkgd[s]==name) return m_systErrors[s];
    }
  }
  return 1;
}

double ProfileLH::getFitParamVal(string name) const {
  if(m_fitSig){
    for(uint s=0; s<m_systNamesFit.size(); ++s){
      if(m_systNamesFit[s]==name) return m_systValues[s];
    }
  }
  else{
    for(uint s=0; s<m_systNamesFitBkgd.size(); ++s){
      if(m_systNamesFitBkgd[s]==name) return m_systValues[s];
    }
  }
  return 0;
}

string ProfileLH::getFitSystName(uint idx){

  if(idx<0) return "";
  if(m_fitSig){
    if(idx>=m_systNamesFit.size()) return ""; 
    return m_systNamesFit[idx];
  }
  else{
    if(idx>=m_systNamesFitBkgd.size()) return ""; 
    return m_systNamesFitBkgd[idx];
  }
  return "";
}

double ProfileLH::getSystSum2(){
  double tot = 0;
  for(uint i=0; i<m_systNames.size(); ++i){
    if(pf_sbd->getFloatFlag(i)) continue;
    tot += pf_sbd->getSystFluctValue(i)*pf_sbd->getSystFluctValue(i);
  }
  return tot;
}

void ProfileLH::fitProfile(){

  if(pf_sbd==NULL){
    printf("ProfileLH::fitProfile, You're trying to fit a profile with a NULL signal/background distribution!!\n"); 
    return;
  }
  
  //Initialize our parameters
  pf_iter = 0; m_status = -1; 
  m_chi2min = -1.0; m_edm = -1.0;
  m_sigScale = m_sigScaleNom; 
  m_sigScaleErr = 0.0;
  pf_sigFloat = m_sigFloat;

  //check to see if the fitting model has changed...
  if((pf_fitSig != m_fitSig) || (m_sigFloat != pf_sigFloat)){
    pf_fitSig = m_fitSig;
    pf_sigFloat = m_sigFloat;
    fillSyst();
    m_init = false;
  }

  //M Owen 16-3-08
  //Only want to vary the systematics in each fit iterations
  //do not want to vary with the stat uncertainties.
  bool histostat = pf_sbd->usingHistoStats();
  pf_sbd->useHistoStats(false);

  if(m_verbose) printf("ProfileLH: %d systematics, %d fit params\n",m_systNames.size(),m_systNamesFit.size());  

  //collect existing parameters  
  for(uint i=0; i<m_systNames.size(); ++i){
    pf_params[i] = pf_sbd->getSystFluctValue(i);
    pf_paramErrs[i] = 0;
  }

  //Record initial conditions...
  m_condI.sig  = pf_sbd->totSignal();
  m_condI.bkgd = pf_sbd->totBkgd();
  m_condI.data = pf_sbd->totData();
  m_condI.chi2 =
    pf_sbd->calculateChi2LLR(pf_fitSig,pf_sigLLR) +
    pf_sbd->calculateSystDiff(pf_params);

  m_condI.chi2_bins = pf_sbd->calculateChi2LLR(pf_fitSig,pf_sigLLR);
  m_condI.chi2_syst = pf_sbd->calculateSystDiff(pf_params);

  char* title = "";
  int ierflg = 0;
  double arglist = -1;
  int narg = 1;
  double step=0.1, bnd1=-4.5, bnd2=4.5;
  //double step=0.2, bnd1=0, bnd2=0;

  if(!m_init || m_sigFloat || m_fitTest){
    MNINIT(0, abs(m_lun), 0);
    
    title = "Collie Profile Likelihood Fitter";
    mnseti_(title,strlen(title));
    MNSETI(title);
    
    ierflg=0;
    char name[250];
    double stval = 0.0;
    for(uint num=1;num<=pf_nParm; ++num){
      sprintf(name,"Par%d",num);
      // bnd1 = 0.0;  // mf 2/10/02 - These two zeros let Minuit used non-bounded
      //  bnd2 = 0.0;  // parameters, at a savings of about 5-10% in time
      MNPARM(num,name,stval,step,bnd1,bnd2,ierflg);
      if(ierflg) {
	fprintf(stderr, "Error on call to mnparm for parameter %s\n",
		name);
	return;
      }
    }
    
    if(m_sigFloat){
      sprintf(name,"SigRate");
      bnd1=1e-5, bnd2=1e5;
      stval = m_sigScaleNom;
      MNPARM(pf_nParm+1,name,stval,step,bnd1,bnd2,ierflg);
      if(ierflg) {
	fprintf(stderr, "Error on call to mnparm for parameter %s\n",
		name);
	return;
      }
    }

    arglist = -1;
    narg = 1;
    title = "SET PRINT";
    MNEXCM(chiFun,title,&arglist,narg,ierflg,NULL);
    
    title = "SET STRAT";
    arglist = 1;  
    MNEXCM(chiFun,title,&arglist,narg,ierflg,NULL);
    
    title = "SET NOW";
    MNEXCM(chiFun,title,&arglist,narg,ierflg,NULL);
    
    title = "SET NOG";
    MNEXCM(chiFun,title,&arglist,narg,ierflg,NULL);
    
    
    arglist = 1;
    title = "SET ERR";
    MNEXCM(chiFun,title,&arglist,narg,ierflg,NULL);
    if(ierflg) {
      fprintf(stderr, "Error on call to mnexcm before fitting\n");
      return;
    }
    m_init = true;      
  }
  
  // Establish the memory of the last Minuit trial point, so that
  // sbd->fillArrays() can operate more efficiently
  pf_lastMinuitTrialPoint.clear();
    
  double arglist2[2];    
  //Protect against failed fits, WF: 8/7/07
  arglist2[0] = m_maxIter;
  arglist2[1] = 1.0;
  narg = 2;
  title = "MINI";
  MNEXCM(chiFun,title,arglist2,narg,ierflg,NULL);
  
  double par3;
  int par4, par5;
  MNSTAT(m_chi2min,m_edm,par3,par4,par5,m_status);

  if(m_status!=3 && m_verbose){
    printf("ProfileLH:  Warning!  MINUIT convergence failing (%d)!\n",m_status);
    //force re-init upon failure
    m_init = false;
  }
  if(m_verbose) printf("ProfileLH: %d MINUIT iterations\n",pf_iter);

  if(pf_iter>0.95*m_maxIter){
    printf("ProfileLH:  MINUIT iteration warning: %d iterations\n",pf_iter);
  }
  
  int ivarbl;
  char chnam[11];
  double val,error;
  m_systErrors.clear();
  m_systValues.clear();

  for(uint num=1;num<=pf_nParm; ++num){
    memset(chnam,0,sizeof(chnam));
    MNPOUT(num,chnam,val,error,bnd1,bnd2,ivarbl);
    m_systErrors.push_back(error);
    m_systValues.push_back(val);
    pf_params[pf_index[num-1]] = val;
    pf_paramErrs[pf_index[num-1]] = error;
  }
  

  if(m_sigFloat){
    memset(chnam,0,sizeof(chnam));
    MNPOUT(pf_nParm+1,chnam,val,error,bnd1,bnd2,ivarbl);
    m_sigScale = val;
    m_sigScaleErr = error;
    pf_sbd->scaleSignal(m_sigScale);    
  }
  
  pf_sbd->fluctuate(pf_params);    


  //Record final conditions...
  m_condF.sig  = pf_sbd->totSignal();
  m_condF.bkgd = pf_sbd->totBkgd();
  m_condF.data = pf_sbd->totData();
  m_condF.chi2 = m_chi2min;
  m_condF.chi2_bins = pf_sbd->calculateChi2LLR(pf_fitSig,pf_sigLLR);
  m_condF.chi2_syst = pf_sbd->calculateSystDiff(pf_params);
  if(m_sigFloat){
    pf_sbd->scaleSignal(1.0/m_sigScale);
    pf_sbd->fluctuate(pf_params);    
  }

  /// Put Histo stats back (M Owen 16-3-08)
  pf_sbd->useHistoStats(histostat);
  //  if(doPrint) print();
  return;
}

// USED TO BE const double* 
std::vector<double>const& ProfileLH::getFitParams() { return pf_params; }

void ProfileLH::print(){
  
  printf("\n****************************************:\n");
  printf("\nProfileLH Fit Results:\n");
  printf("Fit signal: %d, Signal LLR cut: %f, Fit Status: %d\n", m_fitSig, m_sigLLR, m_status);
  printf("MINUIT chi2 min: %f, edm: %f\n",m_chi2min, m_edm);
  printf("Initial Conditions:\n");
  printf("==>Signal: %.3f\n",m_condI.sig);
  printf("==>Bkgd  : %.3f\n",m_condI.bkgd);
  printf("==>Data  : %.0f\n",m_condI.data);
  printf("==>Chi2  : %.3f (%.3f/%.3f)\n",m_condI.chi2,m_condI.chi2_bins,m_condI.chi2_syst);
  printf("Fitted Conditions:\n");
  printf("==>Signal: %.3f\n",m_condF.sig);
  printf("==>Bkgd  : %.3f\n",m_condF.bkgd);
  printf("==>Data  : %.0f\n",m_condF.data);
  printf("==>Chi2  : %.3f (%.3f/%.3f)\n",m_condF.chi2,m_condF.chi2_bins,m_condF.chi2_syst);
  for(uint i=0; i<pf_nParm; ++i){
    printf("==>Fit Param: %14s, Value: %.4f,\t Error: %.4f\n",
	   m_systNames[pf_index[i]].c_str(), 
	   m_systValues[i], 
	   m_systErrors[i]);
  }
  if(m_sigFloat){
    printf("==>Sig Scale: %f, Error: %f\n",m_sigScale,m_sigScaleErr);
  }
  printf("****************************************:\n");
  return;
}

void chiFun(int * npar,double* grad,double *fval,
	    double* par,int * iflag, void (*dummy)()){

  static TrialPoint stat_trialPointDummy;
  if(*iflag==1) {
    /*
     *      Initialize.
     */
    std::cout << " ProfileLH fit function initialize called" << std::endl;
  }
  else if(*iflag==2){
    /*
     *        derivatives...
     */
    /*
    if(pf_fitSig){
      pf_sbd->fluctuate(par);
      for(int i=0; i<pf_nParm; i++) {
	grad[i] = pf_sbd->calculateChi2LLRderivative(i,pf_fitSig,pf_sigLLR,par);
	grad[i] += 2.0*par[i];
      }      
    }
    else{
      for(int i=0; i<pf_nParm; i++) pf_params[pf_index[i]] = par[i];
      pf_sbd->fluctuate(pf_params);
      for(int i=0; i<pf_nParm; i++) {
	grad[i] = pf_sbd->calculateChi2LLRderivative(pf_index[i],pf_fitSig,pf_sigLLR,pf_params);
	grad[i] += 2.0*pf_params[pf_index[i]];
      }
    }
    */
  }    
  else{    
    ++pf_iter;
    
    pf_scl = 1.00;
    if(pf_sigFloat){ 
      pf_scl = par[pf_nParm]; 
      if(pf_scl<1e-5 || std::isnan(pf_scl)!=0 || std::isinf(pf_scl)!=0) pf_scl=1e-5; 
    }

    for(pf_idx=0; pf_idx<pf_nParm; ++pf_idx) pf_params[pf_index[pf_idx]] = par[pf_idx];
    
    pf_sbd->fluctuate(pf_params, pf_lastMinuitTrialPoint);
    if(pf_sigFloat) pf_sbd->scaleSignal(pf_scl);
    *fval = pf_sbd->calculateChi2LLR(pf_fitSig,pf_sigLLR);
    *fval += pf_sbd->calculateSystDiff(pf_params);
    if(pf_sigFloat) pf_sbd->scaleSignal(1.0/pf_scl);    
  } 

  return;
}

void fakeFunctionToSuppressCfortranAndRelatedWarnings() {
  c2fstrv("0","0",0,0);
  f2cstrv("0","0",0,0);
  vkill_trailing("0",0,0,'0');
  num_elem("0",0,0,0); 
  MNUNPT('0');
}
