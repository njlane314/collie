#include "ProfileLH_2D.hh"
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

SigBkgdDist* pf2d_sbd;
std::vector<double> pf2d_params; 
std::vector<uint> pf2d_index;
uint pf2d_nParm;
uint pf2d_iter;
uint pf2d_idx;
bool pf2d_fitSig;
bool pf2d_sigFloat;
double pf2d_sigLLR;
double pf2d_scl1;
double pf2d_scl2;
TrialPoint pf2d_lastMinuitTrialPoint;
#endif // ProfileParams_DEF

///Standard constructor
ProfileLH_2D::ProfileLH_2D(){
  setup();
}

///Constructor with SBD initializer
ProfileLH_2D::ProfileLH_2D(SigBkgdDist* asbd){
  if(asbd==NULL){
    printf("ProfileLH_2D::ProfileLH_2D, the input SBD is NULL!\n");
    return;
  }
  setup();
  setModel(asbd);
}



void ProfileLH_2D::setup(){
  pf2d_sbd       = NULL;
  pf2d_params.clear();  // USED TO BE:  pf2d_params = NULL;
  pf2d_index.clear();
  pf2d_paramErrs.clear();
  pf2d_nParm = 0;
  pf2d_iter = 0;
  pf2d_fitSig   = true;
  pf2d_sigFloat = false;  
  pf2d_sigLLR = 1e6;
  pf2d_scl1 = 1.0;
  pf2d_scl2 = 1.0;

  m_init      = false;
  m_fitTest   = false;
  m_fitSig    = true;
  m_verbose   = false;
  m_sigFloat  = false;  

  m_sigScale1    = 1.0;  
  m_sigScale2    = 1.0;  

  m_sigScaleNom1 = 1.0;
  m_sigScaleNom2 = 1.0;

  m_sigScaleErr2 = 0.0;
  m_sigScaleErr1 = 0.0;

  m_sigLLR      = 1e6;
  m_chi2min     = -1.0;
  m_edm         = -1.0;
  m_lun         = 99;
  m_status      = -1;
  m_maxIter     = 55000;

  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();
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

void ProfileLH_2D::setModel(SigBkgdDist* asbd){  
  pf2d_sbd = asbd;
  pf2d_fitSig = m_fitSig;
  pf2d_sigLLR = m_sigLLR;
  fillSyst();
}

int ProfileLH_2D::getNiterations(){
  return pf2d_iter;
}

void ProfileLH_2D::getMNerrs(){
  if(!m_init) { printf("ProfileLH_2D::getMNerrs, not initialized\n"); return;}
  
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
  
  double errP, errM, eparab, globcc;
  for(uint num=1;num<=pf2d_nParm; ++num){
    arglist[1] = num;
    MNEXCM(chiFun2D,title,arglist,narg,ierflg,NULL);
    MNERRS(num,errP,errM,eparab,globcc);
    m_systErrorsUp.push_back(errP);
    m_systErrorsDn.push_back(errM);
    m_systErrorsParab.push_back(eparab);
    m_systErrorsGlobCC.push_back(globcc);
  }

  return;
}

bool ProfileLH_2D::getEMat(double* emat){
  
  if(!m_init) return false;
  
  //Deriving MINOS errors... 
  double arglist[2];
  arglist[0] = m_maxIter;
  arglist[1] = 1.0;
  char* title = "MINOS";
  int narg = 2;
  int ierflg=0;  
  MNEXCM(chiFun2D,title,arglist,narg,ierflg,NULL);

  double par3;
  int par4, par5;
  int status;
  double chi2min;
  MNSTAT(chi2min,m_edm,par3,par4,par5,status);
  if(status!=3) printf("\nProfileLH2D::getEMat, Status = %d\n\n",status);

  MNEMAT(*emat,pf2d_nParm);  
  
  double errP, errM, eparab, globcc;
  for(uint num=1;num<=pf2d_nParm; ++num){
    arglist[1] = num;
    MNEXCM(chiFun2D,title,arglist,narg,ierflg,NULL);
    MNERRS(num,errP,errM,eparab,globcc);
    m_systErrorsUp.push_back(errP);
    m_systErrorsDn.push_back(errM);
    m_systErrorsParab.push_back(eparab);
    m_systErrorsGlobCC.push_back(globcc);
  }

  if(m_sigFloat){
    arglist[1] = pf2d_nParm+1;
    MNEXCM(chiFun2D,title,arglist,narg,ierflg,NULL);
    MNERRS(pf2d_nParm+1,errP,errM,eparab,globcc);
    
    //    m_sigScaleErr1 = (fabs(errP)+fabs(errM))/2.0;    
    arglist[1] = pf2d_nParm+2;
    MNEXCM(chiFun2D,title,arglist,narg,ierflg,NULL);
    MNERRS(pf2d_nParm+2,errP,errM,eparab,globcc);
    //    m_sigScaleErr2 = (fabs(errP)+fabs(errM))/2.0;    
    //    printf("mnerrs: %f, %f, %f, %f, %f\n",m_sigScaleErr,errP,errM,eparab,globcc);

  }

  return true;
}

void ProfileLH_2D::fillSyst(){
  m_systNames.clear();
  m_systNamesFit.clear();
  m_systNamesFitBkgd.clear();
  m_systNamesFitMap.clear();
  m_systNamesFitBkgdMap.clear();

  if(pf2d_sbd==NULL){
    printf("ProfileLH_2D::fillSyst, You're trying to fit a profile with a NULL signal/background distribution!!\n"); 
    return;
  }
  
  //Make sure all the systematics are properly set
  if(!pf2d_sbd->isInit()){
    pf2d_sbd->varySystematics();    
    pf2d_sbd->setBaselineModel();
  }
  
  vector<string> chanNames = pf2d_sbd->getChannelNames();
  const CollieDistribution* dist;
  string systName ="";
  map<string,int> m_systMap;

  for(vector<string>::iterator its = chanNames.begin(); its!=chanNames.end(); ++its){
    for(uint b = 0; b<pf2d_sbd->getNbkgdDist(*its); ++b){
      
      dist = pf2d_sbd->getBkgdDist(*its,b);
      if(dist!=NULL){
	
	for(int s = 0; s<dist->getNsystematics(); ++s){
	  systName = dist->getSystName(s);
	  if(m_systMap.find(systName)!=m_systMap.end()) continue;
	  else{
	    m_systNames.push_back(systName);	    
	    if(pf2d_sbd->getSystFitFlag(systName)){
	      m_systNamesFitMap[systName] = pf2d_sbd->getSystIndex(systName);
	      m_systNamesFitBkgdMap[systName] = pf2d_sbd->getSystIndex(systName);
	      m_systNamesFit.push_back(systName);
	      m_systNamesFitBkgd.push_back(systName);
	    }
	    else{
	      if(m_verbose) printf("ProfileLH_2D:: Not Fitting %s\n",dist->getSystName(s).c_str());
	    }
	    m_systMap[systName] = pf2d_sbd->getSystIndex(systName);
	  }
	}
      }
    }
  }
  
  for(vector<string>::iterator its = chanNames.begin(); its!=chanNames.end(); ++its){
    for(uint b = 0; b<pf2d_sbd->getNsigDist(*its); ++b){
      
      dist = pf2d_sbd->getSigDist(*its,b);
      if(dist!=NULL){	
	
	for(int s = 0; s<dist->getNsystematics(); ++s){	    
	  systName = dist->getSystName(s);
	  
	  if(m_systMap.find(systName)!=m_systMap.end()) continue;
	  else{
	    m_systNames.push_back(systName);
	    if(pf2d_sbd->getSystFitFlag(systName) && pf2d_fitSig){
	      m_systNamesFitMap[systName] = pf2d_sbd->getSystIndex(systName);
	      m_systNamesFit.push_back(systName);
	    }
	    else{
	      if(m_verbose) printf("ProfileLH_2D:: Not Fitting %s\n",systName.c_str());
	    }
	    m_systMap[systName] = pf2d_sbd->getSystIndex(systName);
	  }
	}
      }
    }
  }

  //our adjustable copy of the syst params
  pf2d_params = std::vector<double>(m_systNames.size(),0.0); 
  pf2d_paramErrs = std::vector<double>(m_systNames.size(),0.0); 
  //  m_systErrors.clear();
  //  m_systValues.clear();  

  //collect existing parameters
  for(unsigned int i=0; i<m_systNames.size(); ++i)  pf2d_params[i] = pf2d_sbd->getSystFluctValue(i);
  
  initFitParams();

  return;
}

void ProfileLH_2D::initFitParams(){
  //Make a storage container with only the values we're fitting...
  //   either background-only or sig+bkgd
  //   for bkgd-only, make index so we know where it fits in master copy.
  //   Allow individual fit parameters to be excluded by the user.
  
  pf2d_index.clear();
  map<string,int>::iterator im;
  
  
  if(m_fitSig){
    pf2d_nParm = m_systNamesFit.size();
    pf2d_index = std::vector<uint>(pf2d_nParm,0); 
    
    //now fill in the syst indices
    for(uint s=0; s<m_systNamesFit.size(); ++s){
      im = m_systNamesFitMap.find(m_systNamesFit[s]);
      pf2d_index[s] = im->second;
    }
  }
  else{
    pf2d_nParm = m_systNamesFitBkgd.size();
    pf2d_index = std::vector<uint>(pf2d_nParm,0); 
    
    //now fill in the syst indices
    for(uint s=0; s<m_systNamesFitBkgd.size(); ++s){
      im = m_systNamesFitBkgdMap.find(m_systNamesFitBkgd[s]);
      pf2d_index[s] = im->second;
    }
  }

  return;
}


bool ProfileLH_2D::getMinosError(uint idx, double& errP, double& errM, double& parab, double& globCC){
  if(m_systErrorsUp.size()==0) getMNerrs();
  if(m_systErrorsUp.size()==0){ printf("ProfileLH_2D::getMinosError, getMinosError failed...\n"); return false;}  //failed..
  if(idx > pf2d_nParm) return false; //failed..
  
  errP = m_systErrorsUp[idx];
  errM = m_systErrorsDn[idx];
  parab = m_systErrorsParab[idx];
  globCC = m_systErrorsGlobCC[idx];
  return true;
}

double ProfileLH_2D::getParamVal(uint p){
  if(p<0 || p>m_systNames.size()) return 0; 
  return pf2d_params[p];
}

double ProfileLH_2D::getParamError(uint p){
  if(p<0 || p>pf2d_paramErrs.size()) return 0; 
  return pf2d_paramErrs[p];
}

int ProfileLH_2D::getFitSystIndex(string name){
  map<string,int>::iterator iter;

  if(name == "Sig1") return pf2d_nParm+1;
  if(name == "Sig2") return pf2d_nParm+2;
  
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

double ProfileLH_2D::getFitParamErr(string name) const {
  
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
  //  printf("ProfileLH2D::getFitParamErr, No such systematic: %s\n",name.c_str());
  return 0;
}

double ProfileLH_2D::getFitParamVal(string name) const {
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
  //  printf("ProfileLH2D::getFitParamValue, No such systematic: %s\n",name.c_str());
  return 0;
}

string ProfileLH_2D::getFitSystName(uint idx){
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

double ProfileLH_2D::getSystSum2(){
  double tot = 0;
  for(uint i=0; i<m_systNames.size(); ++i){
    if(pf2d_sbd->getFloatFlag(i)) continue;
    tot += pf2d_sbd->getSystFluctValue(i)*pf2d_sbd->getSystFluctValue(i);
  }
  return tot;
}

void ProfileLH_2D::fitProfile(){

  if(pf2d_sbd==NULL){
    printf("ProfileLH_2D::fitProfile, You're trying to fit a profile with a NULL signal/background distribution!!\n"); 
    return;
  }
  
  //Initialize our parameters
  pf2d_iter = 0; m_status = -1; 
  m_chi2min = -1.0; m_edm = -1.0;
  m_sigScale1 = m_sigScaleNom1; 
  m_sigScaleErr1 = 0.0;

  m_sigScale2 = m_sigScaleNom2; 
  m_sigScaleErr2 = 0.0;

  pf2d_sigFloat = m_sigFloat;

  //check to see if the fitting model has changed...
  if((pf2d_fitSig != m_fitSig) || (m_sigFloat != pf2d_sigFloat)){
    pf2d_fitSig = m_fitSig;
    pf2d_sigFloat = m_sigFloat;
    fillSyst();
    m_init = false;
  }

  //M Owen 16-3-08
  //Only want to vary the systematics in each fit iterations
  //do not want to vary with the stat uncertainties.
  bool histostat = pf2d_sbd->usingHistoStats();
  pf2d_sbd->useHistoStats(false);

  if(m_verbose) printf("ProfileLH_2D: %d systematics, %d fit params\n",m_systNames.size(),m_systNamesFit.size());  

  //collect existing parameters  
  for(uint i=0; i<m_systNames.size(); ++i){
    pf2d_params[i] = pf2d_sbd->getSystFluctValue(i);
    pf2d_paramErrs[i] = 0;
  }

  //Record initial conditions...
  m_condI.sig  = pf2d_sbd->totSignal();
  m_condI.bkgd = pf2d_sbd->totBkgd();
  m_condI.data = pf2d_sbd->totData();
  m_condI.chi2 =
    pf2d_sbd->calculateChi2LLR(pf2d_fitSig,pf2d_sigLLR) +
    pf2d_sbd->calculateSystDiff(pf2d_params);

  m_condI.chi2_bins = pf2d_sbd->calculateChi2LLR(pf2d_fitSig,pf2d_sigLLR);
  m_condI.chi2_syst = pf2d_sbd->calculateSystDiff(pf2d_params);


  char* title = "";
  int ierflg = 0;
  double arglist = -1;
  int narg = 1;
  double step=0.1, bnd1=-4.5, bnd2=4.5;

  if(!m_init || m_sigFloat || m_fitTest){
    MNINIT(0, abs(m_lun), 0);
    
    title = "Collie Profile Likelihood Fitter";
    mnseti_(title,strlen(title));
    MNSETI(title);
    
    ierflg=0;
    char name[250];
    double stval = 0.0;
    for(uint num=1;num<=pf2d_nParm; ++num){
      sprintf(name,"Par%d",num);
      //      bnd1 = 0.0;  // mf 2/10/02 - These two zeros let Minuit used non-bounded
      //      bnd2 = 0.0;  // parameters, at a savings of about 5-10% in time
      MNPARM(num,name,stval,step,bnd1,bnd2,ierflg);
      if(ierflg) {
	fprintf(stderr, "Error on call to mnparm for parameter %s\n",
		name);
	return;
      }
    }
    
    if(m_sigFloat){
      sprintf(name,"SigRate1");
      bnd1=1e-5, bnd2=1e5;
      stval = m_sigScaleNom1;
      MNPARM(pf2d_nParm+1,name,stval,step,bnd1,bnd2,ierflg);

      sprintf(name,"SigRate2");
      bnd1=1e-5, bnd2=1e5;
      stval = m_sigScaleNom2;
      MNPARM(pf2d_nParm+2,name,stval,step,bnd1,bnd2,ierflg);
      if(ierflg) {
	fprintf(stderr, "Error on call to mnparm for parameter %s\n",
		name);
	return;
      }
    }

    arglist = -1;
    narg = 1;
    title = "SET PRINT";
    MNEXCM(chiFun2D,title,&arglist,narg,ierflg,NULL);
    
    title = "SET STRAT";
    arglist = 1;  
    MNEXCM(chiFun2D,title,&arglist,narg,ierflg,NULL);
    
    title = "SET NOW";
    MNEXCM(chiFun2D,title,&arglist,narg,ierflg,NULL);
    
    title = "SET NOG";
    MNEXCM(chiFun2D,title,&arglist,narg,ierflg,NULL);
    
    
    arglist = 1;
    title = "SET ERR";
    MNEXCM(chiFun2D,title,&arglist,narg,ierflg,NULL);
    if(ierflg) {
      fprintf(stderr, "Error on call to mnexcm before fitting\n");
      return;
    }
    m_init = true;      
  }
  
  // Establish the memory of the last Minuit trial point, so that
  // sbd->fillArrays() can operate more efficiently
  pf2d_lastMinuitTrialPoint.clear();
  
  
  double arglist2[2];  
  //Protect against failed fits, WF: 8/7/07
  arglist2[0] = m_maxIter;
  arglist2[1] = 1.0;
  narg = 2;
  title = "MINI";
  MNEXCM(chiFun2D,title,arglist2,narg,ierflg,NULL);

  double par3;
  int par4, par5;
  MNSTAT(m_chi2min,m_edm,par3,par4,par5,m_status);

  if(m_status!=3 && m_verbose){
    printf("ProfileLH_2D:  Warning!  MINUIT convergence failing (%d)!\n",m_status);
    //force re-init upon failure
    m_init = false;
  }
  if(m_verbose) printf("ProfileLH_2D: %d MINUIT iterations\n",pf2d_iter);

  if(pf2d_iter>0.95*m_maxIter){
    printf("ProfileLH_2D:  MINUIT iteration warning: %d iterations\n",pf2d_iter);
  }
  
  int ivarbl;
  char chnam[11];
  double val,error;
  m_systErrors.clear();
  m_systValues.clear();

  for(uint num=1;num<=pf2d_nParm; ++num){
    memset(chnam,0,sizeof(chnam));
    MNPOUT(num,chnam,val,error,bnd1,bnd2,ivarbl);
    m_systErrors.push_back(error);
    m_systValues.push_back(val);
    pf2d_params[pf2d_index[num-1]] = val;
    pf2d_paramErrs[pf2d_index[num-1]] = error;
  }
  

  if(m_sigFloat){
    memset(chnam,0,sizeof(chnam));
    MNPOUT(pf2d_nParm+1,chnam,val,error,bnd1,bnd2,ivarbl);
    m_sigScale1 = val;
    m_sigScaleErr1 = error;
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0),m_sigScale1, pf2d_params);    

    memset(chnam,0,sizeof(chnam));
    MNPOUT(pf2d_nParm+2,chnam,val,error,bnd1,bnd2,ivarbl);
    m_sigScale2 = val;
    m_sigScaleErr2 = error;
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1),m_sigScale2, pf2d_params);    
  }
  
  pf2d_sbd->fluctuate(pf2d_params);    


  //Record final conditions...
  m_condF.sig  = pf2d_sbd->totSignal();
  m_condF.bkgd = pf2d_sbd->totBkgd();
  m_condF.data = pf2d_sbd->totData();
  m_condF.chi2 = m_chi2min;
  m_condF.chi2_bins = pf2d_sbd->calculateChi2LLR(pf2d_fitSig,pf2d_sigLLR);
  m_condF.chi2_syst = pf2d_sbd->calculateSystDiff(pf2d_params);
  if(m_sigFloat){
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0),1.0/m_sigScale1, pf2d_params);
    pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1),1.0/m_sigScale2, pf2d_params);
    pf2d_sbd->fluctuate(pf2d_params);    
  }

  /// Put Histo stats back (M Owen 16-3-08)
  pf2d_sbd->useHistoStats(histostat);
  //  if(doPrint) print();
  return;
}

// USED TO BE const double* 
std::vector<double>const& ProfileLH_2D::getFitParams() { return pf2d_params; }

void ProfileLH_2D::print(){
  
  printf("\n****************************************:\n");
  printf("\nProfileLH_2D Fit Results:\n");
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
  for(uint i=0; i<pf2d_nParm; ++i){
    printf("==>Fit Param: %14s, Value: %.4f,\t Error: %.4f\n",
	   m_systNamesFit[pf2d_index[i]].c_str(), 
	   m_systValues[i], 
	   m_systErrors[i]);
  }
  if(m_sigFloat){
    printf("==>Sig1 Scale: %f, Error: %f\n",m_sigScale1,m_sigScaleErr1);
    printf("==>Sig2 Scale: %f, Error: %f\n",m_sigScale2,m_sigScaleErr2);
  }
  printf("****************************************:\n");
  return;
}

void chiFun2D(int * npar,double* grad,double *fval,
	      double* par,int * iflag, void (*dummy)()){

  //  static TrialPoint stat_trialPointDummy;
  if(*iflag==1) {
    /*
     *      Initialize.
     */
    std::cout << " ProfileLH_2D fit function initialize called" << std::endl;
  }
  else if(*iflag==2){
    /*
     *        derivatives...
     */
    /*
    if(pf2d_fitSig){
      pf2d_sbd->fluctuate(par);
      for(int i=0; i<pf2d_nParm; ++i) {
	grad[i] = pf2d_sbd->calculateChi2LLRderivative(i,pf2d_fitSig,pf2d_sigLLR,par);
	grad[i] += 2.0*par[i];
      }      
    }
    else{
      for(int i=0; i<pf2d_nParm; ++i) pf2d_params[pf2d_index[i]] = par[i];
      pf2d_sbd->fluctuate(pf2d_params);
      for(int i=0; i<pf2d_nParm; ++i) {
	grad[i] = pf2d_sbd->calculateChi2LLRderivative(pf2d_index[i],pf2d_fitSig,pf2d_sigLLR,pf2d_params);
	grad[i] += 2.0*pf2d_params[pf2d_index[i]];
      }
    }
    */
  }    
  else{    
    ++pf2d_iter;
    
    pf2d_scl1 = 1.00;
    pf2d_scl2 = 1.00;
    if(pf2d_sigFloat){ 

      pf2d_scl1 = par[pf2d_nParm]; 
      if(pf2d_scl1<1e-5 || std::isnan(pf2d_scl1)!=0 || std::isinf(pf2d_scl1)!=0) pf2d_scl1=1e-5; 

      pf2d_scl2 = par[pf2d_nParm+1]; 
      if(pf2d_scl2<1e-5 || std::isnan(pf2d_scl2)!=0 || std::isinf(pf2d_scl2)!=0) pf2d_scl2=1e-5; 

    }

    for(pf2d_idx=0; pf2d_idx<pf2d_nParm; ++pf2d_idx) pf2d_params[pf2d_index[pf2d_idx]] = par[pf2d_idx];
    
    pf2d_sbd->fluctuate(pf2d_params, pf2d_lastMinuitTrialPoint);

    if(pf2d_sigFloat){
      pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0),pf2d_scl1,pf2d_params);
      pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1),pf2d_scl2,pf2d_params);
    }

    *fval = pf2d_sbd->calculateChi2LLR(pf2d_fitSig,pf2d_sigLLR);
    *fval += pf2d_sbd->calculateSystDiff(pf2d_params);

    if(pf2d_sigFloat){
      pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(0),1.0/pf2d_scl1,pf2d_params);
      pf2d_sbd->scaleSignal(pf2d_sbd->getSigName(1),1.0/pf2d_scl2,pf2d_params);
    }
  } 

  return;
}

void fakeFunctionToSuppressCfortranAndRelatedWarnings2D() {
  c2fstrv("0","0",0,0);
  f2cstrv("0","0",0,0);
  vkill_trailing("0",0,0,'0');
  num_elem("0",0,0,0); 
  MNUNPT('0');
}
