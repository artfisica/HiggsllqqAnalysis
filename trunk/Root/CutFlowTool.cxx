#define __CutFlowTool_C__

#include "HiggsllqqAnalysis/CutFlowTool.h"

ClassImp(Analysis::CutFlowTool);

void Analysis::CutFlowTool::addCut(TString cut_name)
{
   UInt_t index = cuts_map_.size();
   cuts_map_[index] = cut_name;
   cut_npass_.push_back(0.);
   return;
}

void Analysis::CutFlowTool::addCutCounter(Int_t cut_index, Double_t nevents)
{
   for (int i = 0; i < cut_index + 1; i++) {
      double old_counter = cut_npass_.at(i);
      cut_npass_.at(i) = old_counter + nevents;
   }
   return;
}

void Analysis::CutFlowTool::resetCounter()
{
   std::map<UInt_t, TString>::iterator it;
   for (it = cuts_map_.begin(); it != cuts_map_.end(); it++) {
      UInt_t index = (*it).first;
      cut_npass_.at(index) = 0.;
   }
   return;
}

void Analysis::CutFlowTool::print()
{
   Info("CutFlowTool", "Printing %s cutflow:", name_.Data());
   Info("CutFlowTool", "|	CUT		| EVENTS\t| ABS. EFF.\t| REL. EFF. |");
   std::map<UInt_t, TString>::iterator it;
   for (it = cuts_map_.begin(); it != cuts_map_.end(); it++) {
      int index = (*it).first;
      int index_minus_one = (index > 0) ? (index - 1) : index;
      if (cut_npass_.at(index_minus_one) != 0 && cut_npass_.at(0) != 0) {
         Info("CutFlowTool", "| %15s | %lf\t| %lf\t | %lf |", ((*it).second).Data(), cut_npass_.at(index), (double)cut_npass_.at(index) / (double)cut_npass_.at(0), (double)cut_npass_.at(index) / (double)cut_npass_.at(index_minus_one));
      } else {
         std::cout << (*it).second << ": " << cut_npass_.at(index) << ", efficiencies could not be computed since the denominators are 0" << std::endl;
      }
   }

   return;
}
