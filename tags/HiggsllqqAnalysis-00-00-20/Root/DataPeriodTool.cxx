#define __DataPeriodTool_C__

#include "HiggsllqqAnalysis/DataPeriodTool.h"

std::pair<UInt_t, UInt_t> DataPeriodTool::GetPeriodRange(TString period)
{
   std::pair<UInt_t, UInt_t> result;

   if (period == "y2011_A") {
      result = std::make_pair(177531, 177965);
   } else if (period == "y2011_B") {
      result = std::make_pair(177986, 178109);
   } else if (period == "y2011_D") {
      result = std::make_pair(179710, 180481);
   } else if (period == "y2011_E") {
      result = std::make_pair(180614, 180776);
   } else if (period == "y2011_F") {
      result = std::make_pair(182013, 182519);
   } else if (period == "y2011_G") {
      result = std::make_pair(182726, 183462);
   } else if (period == "y2011_H") {
      result = std::make_pair(183544, 184169);
   } else if (period == "y2011_I") {
      result = std::make_pair(185353, 186493);
   } else if (period == "y2011_J") {
      result = std::make_pair(186516, 186755);
   } else if (period == "y2011_K") {
      result = std::make_pair(186873, 187815);
   } else if (period == "y2011_L") {
      result = std::make_pair(188902, 190343);
   } else if (period == "y2011_M") {
      result = std::make_pair(190503, 191933);
   } else if (period == "y2012_A") {
      result = std::make_pair(200804, 201556);
   } else if (period == "y2012_B") {
      result = std::make_pair(202660, 205113); // update it often!!!
   }

   return result;
}

TString DataPeriodTool::GetPeriodName(UInt_t run)
{
   TString dummy = "y2011_none";

   // list of periods (MUST BE ORDERED!!!)
   std::vector<TString> periods;
   periods.push_back("y2011_A");
   periods.push_back("y2011_B");
   periods.push_back("y2011_C");
   periods.push_back("y2011_D");
   periods.push_back("y2011_E");
   periods.push_back("y2011_F");
   periods.push_back("y2011_G");
   periods.push_back("y2011_H");
   periods.push_back("y2011_I");
   periods.push_back("y2011_J");
   periods.push_back("y2011_K");
   periods.push_back("y2011_L");
   periods.push_back("y2011_M");
   periods.push_back("y2012_A");
   periods.push_back("y2012_B");

   for (UInt_t i = 0; i < periods.size(); i++) {
      TString this_period = periods.at(i);

      if (run <= DataPeriodTool::GetPeriodRange(this_period).second) {
         return this_period;
      }
   }

   return dummy;
}
