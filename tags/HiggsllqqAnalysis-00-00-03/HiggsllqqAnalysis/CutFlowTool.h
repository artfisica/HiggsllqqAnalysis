#ifndef __CutFlowTool_h__
#define __CutFlowTool_h__


/** @class CutFlowTool CutFlowTool.h

    Functions to deal with cut flows.

    @author Camilla Maiani <camilla.maiani@cern.ch>
    @author Valerio Ippolito <valerio.ippolito@cern.ch>
    @date 02/11/2011
*/


#include <iostream>
#include <vector>
#include <map>
#include <TString.h>
#include <TObject.h>

namespace Analysis {

   class CutFlowTool : public TObject {
   public:

      CutFlowTool(TString name = "") {
         name_ = name;
         cuts_map_.clear();
         cut_npass_.clear();
      };
      ~CutFlowTool() {};

      TString getName() {
         return name_;
      }
      std::map<UInt_t, TString> getCutFlowToolMap() {
         return cuts_map_;
      }
      std::vector<Double_t> getNPassVec() {
         return cut_npass_;
      }
      Double_t getNPass(UInt_t index) {
         return cut_npass_[index];
      }
      TString getCutName(UInt_t index) {
         return cuts_map_[index];
      }

      void addCut(TString cut_name);
      void addCutCounter(Int_t cut_index, Double_t nevents);
      void resetCounter();
      void print();

   private:
      TString name_;
      std::map<UInt_t, TString> cuts_map_;
      std::vector<Double_t> cut_npass_;

      ClassDef(CutFlowTool, 1);
   };
};

#endif
