#include "HiggsllqqAnalysis/HiggsllqqAnalysis.h"
#include <Cintex/Cintex.h>
#include <TTreeCache.h>

void print_usage(char *app_name)
{
  std::cout << "Usage:" << std::endl;
  std::cout << app_name << " [--MUID] [--GSF] [--JET] [--noSmearing] [--useTopoIso] [--analysis <rel_17|rel_17_2>] [--input <filename.txt>] [--output <filename.root>] [--MV1c]" << std::endl;
}

int main(int argc, char **argv)
{
  ROOT::Cintex::Cintex::Enable();
  
  // option parser
  Bool_t doSTACO(kTRUE);
  Bool_t doGSF(kFALSE);
  Bool_t doJET(kFALSE);
  Bool_t doMV1c(kFALSE);
  Bool_t doSmearing(kTRUE);
  Bool_t useTopoIso(kFALSE);
  TString output_filename("output_test.root");
  TString analysis_version("rel_17");
  TString input_filename("input.txt");
  
  Bool_t error_in_options(kFALSE);
  
  for (int i = 1; i < argc; i++)
    {
      if (TString(argv[i]) == "--GSF")             doGSF      = kTRUE;
      else if (TString(argv[i]) == "--MUID")       doSTACO    = kFALSE;
      else if (TString(argv[i]) == "--JET")        doJET      = kTRUE;
      else if (TString(argv[i]) == "--MV1c")       doMV1c     = kTRUE;
      else if (TString(argv[i]) == "--noSmearing") doSmearing = kFALSE;
      else if (TString(argv[i]) == "--useTopoIso") useTopoIso = kTRUE;
      else if (TString(argv[i]) == "--analysis") {
	if (i + 1 < argc)
	  {
	    analysis_version = argv[++i];
	  } else {
	  error_in_options = kTRUE;
	}
      }
      else if (TString(argv[i]) == "--input")
	{
	  if (i + 1 < argc)
	    {
	      input_filename = argv[++i];
	    } 
	  else
	    {
	      error_in_options = kTRUE;
	    }
	} else if (TString(argv[i]) == "--output")
	{
	  if (i + 1 < argc)
	    {
	      output_filename = argv[++i];
	    } 
	  else
	    {
	      error_in_options = kTRUE;
	    }
	} 
      else
	{
	  error_in_options = kTRUE;
	}
    }
   
  if (error_in_options)
    {
      print_usage(argv[0]);
      return 1;
    } 
  else
    {
      std::cout << argv[0] << " called with options:"           << std::endl;
      std::cout << "  analysis_version = " << analysis_version  << std::endl;
      std::cout << "           doSTACO = " << doSTACO           << std::endl;
      std::cout << "             doGSF = " << doGSF             << std::endl;
      std::cout << "             doJET = " << doJET             << std::endl;
      std::cout << "            doMV1c = " << doMV1c            << std::endl;
      std::cout << "        doSmearing = " << doSmearing        << std::endl;
      std::cout << "        useTopoIso = " << useTopoIso        << std::endl;
      std::cout << "    input_filename = " << input_filename    << std::endl;
      std::cout << "   output_filename = " << output_filename   << std::endl;
      std::cout << std::endl;
    }
  ///
  
  TChain *theChain = new TChain("physics");
  
  std::ifstream input;
  input.open(input_filename.Data());
  
  std::string ntuple;
  
  input >> ntuple;
  
  while (!input.eof())
    {
      theChain->Add(ntuple.c_str());
      input >> ntuple;
    }
  
  // TTreeCache
  //TTreeCache::SetLearnEntries(100);
  theChain->SetCacheSize(10000000); // 10 Mb
  
  HiggsllqqAnalysis *analysis = new HiggsllqqAnalysis();
  analysis->setAnalysisVersion(analysis_version);
  analysis->setMuonFamily((doSTACO)   ? Muon::STACO : Muon::MUID);
  analysis->setElectronFamily((doGSF) ? Electron::GSF : Electron::noGSF);
  analysis->setJetFamily((doJET)      ? 1 : 0);
  analysis->setJetbTagger((doMV1c)    ? 1 : 0);
  analysis->setSmearing(doSmearing);
  analysis->setTopoIso(useTopoIso);
  analysis->setOutputFile(output_filename);
  
  theChain->Process(analysis);
  
  theChain->PrintCacheStats();
  
  delete analysis;
  delete theChain;
  
  return 0;
}
