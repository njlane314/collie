#include <TFile.h>
#include <TTree.h>
#include <CrossSectionLimit.hh>
#include <CollieLoader.hh>
#include <FitTest.hh>
#include <CLfast.hh>
#include <CLsyst.hh>
#include <CLfit.hh>
#include <CLfit2.hh>
#include <sys/time.h>


void calcLimit(char* outFile, char* inList, char* m){

  int mass = atoi(m);

  printf("\n\ncollieLimitCalc.exe\n");
  printf("Input List: %s\n",inList);
  printf("Output File: %s\n",outFile);
  printf("Mass Point: %d\n",mass);

  CollieLoader loaders[300];
  string chanNames[300];
  string fileNames[300];
  int nld = 0;
  ifstream streamIn(inList);
  if(!streamIn){    
    cout << "Error: Could not open " << inList << endl;
    return;
  } 
  
  bool ok = true;
  char fname[1024];  
  char options[1024];
  while(!streamIn.eof()){
    if(!(streamIn >> fname)) continue;
    cout << "Reading: " << fname << endl;
    
    TFile* ftest = new TFile(fname);
    TList* aList = ftest->GetListOfKeys();
    if(aList->GetEntries()!=1){
      printf("Incorrect key length for %s\n",fname);
      return;
    }
    chanNames[nld] = aList->At(0)->GetName();
    TString name(fname);
    fileNames[nld] = name.Data();
    aList->Delete();
    ftest->Close();
    ftest->Delete();
    
    
    sprintf(options,"name='%s'",chanNames[nld].c_str());
    if (!loaders[nld].open(fname,options)) {
      std::cout << "Failed to open " << fname << " using " << options << "!\n";
      ok = false;
    }
    
    nld++;
  }
  
  printf("\n************************************************\n");
  printf("Collie Example Limit Calculation\n");
  printf("%d channel(s) available\n",nld);
  printf("************************************************\n");

  if(!ok) return;

