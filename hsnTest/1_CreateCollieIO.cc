/* This script takes in a root file containing histograms for the signal, background and data distributions and creates a CollieIOfile that can be read by the successive scripts.

The input file may contain several histograms for signal, for each of the combination of parameters used in the simulation (e.g. mass and mixing angle).

For each combination of parameters a createMassPoint method is called, which saves a mass point in the output file. Each mass point takes in a signal, data and vector of background histograms.
*/

#include "CollieIOFile.hh"
#include "TRandom.h"
#include "TClass.h"
#include <string>
#include <sstream>

int main(int argc, char* argv[])
{
  std::string label = std::string(argv[1]);
  int massName = atoi(argv[2]);
  int run = atoi(argv[3]);

  stringstream inString;
  inString << "/pc2014-data5/nlane/hnl/collie/hsnTest/Storage/InputFiles/" << label << "_m" << massName  << "_r" << run << ".root";
  
  stringstream outString;
  outString << "/pc2014-data5/nlane/hnl/collie/hsnTest/Storage/CollieFiles/CollieIOfile_" << label << "_m" << massName << "_r" << run << ".root";
  
  stringstream chanName;
  chanName  << "HNLChannel_r" << run;

  // Open root input file
  TFile infile(inString.str().c_str());
  TH1D* data = (TH1D*)infile.Get("data");
  TH1D* bkgd1 = (TH1D*)infile.Get("bkg_EXT");
  TH1D* bkgd2 = (TH1D*)infile.Get("bkg_overlay");
  TH1D* bkgd3 = (TH1D*)infile.Get("bkg_dirt");

  // Specify outputfile and channel name 
  CollieIOFile* cfile = new CollieIOFile();
  cfile->initFile(outString.str().c_str(), chanName.str().c_str());  
  cfile->setNoviceFlag(false);

  // Define your input histograms
  double Xmin = 0; 
  double Xmax = 1;
  int Nbins = 5; 
  cfile->setInputHist(Xmin,Xmax,Nbins);
  cfile->setRebin(1);

  // Define background names
  vector<string> bkgdNames;
  bkgdNames.push_back("ext3");
  bkgdNames.push_back("nu3");
  bkgdNames.push_back("dirt3");
  cfile->createChannel(bkgdNames);
  
  // Vector of backgrounds
  vector<TH1D*> vbkgd;
  vbkgd.push_back(bkgd1);
  vbkgd.push_back(bkgd2);
  vbkgd.push_back(bkgd3);
  // Vector of useless alpha parameters
  // (for smoothing, which we are not using)
  vector<double> valpha;

  // Create a parameter point to analyze
  int mass[1] = {massName};
  for (int m : mass)
  {
    std::string histName = std::string("signal");
    TString tHistName(histName);
    TH1D* signal = (TH1D*)infile.Get(tHistName);
    cfile->createMassPoint(m, data, signal, -1, vbkgd, valpha);
    cfile->createFlatSigSystematic("Flux",0.3,0.3,m);
    cfile->createFlatSigSystematic("SigDect",0.15,0.15,m);
    cfile->createFlatSigSystematic("POT",0.02,0.02,m);

    TH1D* ppfx_up = (TH1D*)infile.Get("ppfx_uncertainty_frac");
    TH1D* ppfx_down = (TH1D*)infile.Get("ppfx_uncertainty_frac");
    TH1D* gen_up = (TH1D*)infile.Get("Genie_uncertainty_frac");
    TH1D* gen_down = (TH1D*)infile.Get("Genie_uncertainty_frac");
    std::cout<<"runnging0"<<std::endl;
    cfile->createFlatBkgdSystematic(0,"POT",0.02,0.02,m);
    cfile->createFlatBkgdSystematic(1,"POT",0.02,0.02,m);
    cfile->createFlatBkgdSystematic(2,"POT",0.02,0.02,m);
    //cfile->createFlatBkgdSystematic(0,"bkg",0.2,0.2,m);
    //cfile->createFlatBkgdSystematic(1,"xsec",0.2,0.2,m);
    cfile->createBkgdSystematic(1, "genie1", gen_up, gen_down, m);
    //cfile->setLogNormalFlag("genie",true,m);
    cfile->createBkgdSystematic(1, "ppfx1", ppfx_up, ppfx_down, m);
    // cfile->setLogNormalFlag("ppfx",true,m);
    cfile->createFlatBkgdSystematic(1,"dectsys1",0.5,0.5,m);
    // cfile->setLogNormalFlag("dectsys",true,m);

    cfile->createFlatBkgdSystematic(2,"dirtnorm",1,1,m);
    cfile->setBkgdFloatFlag(2,"dirtnorm",true,m);
  }
    
  cfile->storeFile();
}

