#ifndef __DataPeriodTool_h__
#define __DataPeriodTool_h__

/** @class DataPeriodTool DataPeriodTool.h

    Functions to deal with data periods and run ranges.

    @author Valerio Ippolito <valerio.ippolito@cern.ch>
    @date 03/11/2011
*/

#include <string>
#include <iostream>
#include <utility>
#include <TString.h>

namespace DataPeriodTool {
   /// returns a std::pair where first and second are the first and last run number within the given period
   std::pair <UInt_t, UInt_t> GetPeriodRange(TString period);
   TString GetPeriodName(UInt_t run);
}

#endif
