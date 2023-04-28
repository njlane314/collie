/* This script takes in a root file containing histograms for the signal, background and data distributions and creates a CollieIOfile that can be read by the successive scripts.

The input file may contain several histograms for signal, for each of the combination of parameters used in the simulation (e.g. mass and mixing angle).

For each combination of parameters a createMassPoint method is called, which saves a mass point in the output file. Each mass point takes in a signal, data and vector of background histograms.
*/

#include "CollieIOFile.hh"
#include "TRandom.h"

int main()
{
  // Open root input file
  TFile infile("../data/e_like/HNL_e_450_all_errors_5e-8.root");
  TH1D* data = (TH1D*)infile.Get("Data");
  TH1D* bkgd1 = (TH1D*)infile.Get("Background");
  TH1D* signal = (TH1D*)infile.Get("Signal");

  // Specify outputfile and channel name 
  CollieIOFile* cfile = new CollieIOFile();
  cfile->initFile("output/sbnd/bdt_CollieIOfile.root", "sbnd_450MeV");  

  // Define your input histograms
  
  double Xmin = 0.0; 
  //double Xmax = 3.2;
  double Xmax = 1.0;
  //int Nbins = 100;
  int Nbins = 20;
  cfile->setInputHist(Xmin,Xmax,Nbins);
  cfile->setRebin(1);

  // Define background names
  vector<string> bkgdNames;
  bkgdNames.push_back("Background");
  cfile->createChannel(bkgdNames);
  
  // Vector of backgrounds
  vector<TH1D*> vbkgd;
  vbkgd.push_back(bkgd1);
  // Vector of useless alpha parameters
  // (for smoothing, which we are not using)
  vector<double> valpha;

  cfile->createMassPoint(0.450, data, signal, -1, vbkgd, valpha);
  // Create a parameter point to analyze
  /*
  int mass[1] = {0.310};
  for (int m : mass)
  {
    std::string histName = std::string("Invmass_signal");
    TString tHistName(histName);
    TH1D* signal = (TH1D*)infile.Get(tHistName);
    cfile->createMassPoint(m, data, signal, -1, vbkgd, valpha);
    std::cout << "Run through mass loop" << std::endl;
  }
  */  
  cfile->storeFile();
}
