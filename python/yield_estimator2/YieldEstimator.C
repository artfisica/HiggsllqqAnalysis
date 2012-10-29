#define YieldEstimator_cxx
// The class definition in YieldEstimator.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("YieldEstimator.C")
// Root > T->Process("YieldEstimator.C","some options")
// Root > T->Process("YieldEstimator.C+")
//

#include "YieldEstimator.h"
#include <TH2.h>
#include <TStyle.h>



Bool_t YieldEstimator::fillMCMap(TString release)
{
   if (release == "mc12a") {
      m_sample[160152] = YieldSample("PowhegPythia8_AU2CT10_ggH110_ZZ4lep", kTRUE, 110, MC::ggF);
      m_sample[160153] = YieldSample("PowhegPythia8_AU2CT10_ggH115_ZZ4lep", kTRUE, 115, MC::ggF);
      m_sample[160154] = YieldSample("PowhegPythia8_AU2CT10_ggH120_ZZ4lep", kTRUE, 120, MC::ggF);
      m_sample[160155] = YieldSample("PowhegPythia8_AU2CT10_ggH125_ZZ4lep", kTRUE, 125, MC::ggF);
      m_sample[160156] = YieldSample("PowhegPythia8_AU2CT10_ggH130_ZZ4lep", kTRUE, 130, MC::ggF);
      m_sample[160157] = YieldSample("PowhegPythia8_AU2CT10_ggH135_ZZ4lep", kTRUE, 135, MC::ggF);
      m_sample[160158] = YieldSample("PowhegPythia8_AU2CT10_ggH140_ZZ4lep", kTRUE, 140, MC::ggF);
      m_sample[160159] = YieldSample("PowhegPythia8_AU2CT10_ggH145_ZZ4lep", kTRUE, 145, MC::ggF);
      m_sample[160160] = YieldSample("PowhegPythia8_AU2CT10_ggH150_ZZ4lep", kTRUE, 150, MC::ggF);
      m_sample[160161] = YieldSample("PowhegPythia8_AU2CT10_ggH155_ZZ4lep", kTRUE, 155, MC::ggF);
      m_sample[160162] = YieldSample("PowhegPythia8_AU2CT10_ggH160_ZZ4lep", kTRUE, 160, MC::ggF);
      m_sample[160163] = YieldSample("PowhegPythia8_AU2CT10_ggH165_ZZ4lep", kTRUE, 165, MC::ggF);
      m_sample[160164] = YieldSample("PowhegPythia8_AU2CT10_ggH170_ZZ4lep", kTRUE, 170, MC::ggF);
      m_sample[160165] = YieldSample("PowhegPythia8_AU2CT10_ggH175_ZZ4lep", kTRUE, 175, MC::ggF);
      m_sample[160166] = YieldSample("PowhegPythia8_AU2CT10_ggH180_ZZ4lep", kTRUE, 180, MC::ggF);
      m_sample[160167] = YieldSample("PowhegPythia8_AU2CT10_ggH185_ZZ4lep", kTRUE, 185, MC::ggF);
      m_sample[160168] = YieldSample("PowhegPythia8_AU2CT10_ggH190_ZZ4lep", kTRUE, 190, MC::ggF);
      m_sample[160169] = YieldSample("PowhegPythia8_AU2CT10_ggH195_ZZ4lep", kTRUE, 195, MC::ggF);
      m_sample[160170] = YieldSample("PowhegPythia8_AU2CT10_ggH200_ZZ4lep", kTRUE, 200, MC::ggF);
      m_sample[160171] = YieldSample("PowhegPythia8_AU2CT10_ggH220_ZZ4lep", kTRUE, 220, MC::ggF);
      m_sample[160172] = YieldSample("PowhegPythia8_AU2CT10_ggH240_ZZ4lep", kTRUE, 240, MC::ggF);
      m_sample[160173] = YieldSample("PowhegPythia8_AU2CT10_ggH260_ZZ4lep", kTRUE, 260, MC::ggF);
      m_sample[160174] = YieldSample("PowhegPythia8_AU2CT10_ggH280_ZZ4lep", kTRUE, 280, MC::ggF);
      m_sample[160175] = YieldSample("PowhegPythia8_AU2CT10_ggH300_ZZ4lep", kTRUE, 300, MC::ggF);
      m_sample[160176] = YieldSample("PowhegPythia8_AU2CT10_ggH320_ZZ4lep", kTRUE, 320, MC::ggF);
      m_sample[160177] = YieldSample("PowhegPythia8_AU2CT10_ggH340_ZZ4lep", kTRUE, 340, MC::ggF);
      m_sample[160178] = YieldSample("PowhegPythia8_AU2CT10_ggH360_ZZ4lep", kTRUE, 360, MC::ggF);
      m_sample[160179] = YieldSample("PowhegPythia8_AU2CT10_ggH380_ZZ4lep", kTRUE, 380, MC::ggF);
      m_sample[160180] = YieldSample("PowhegPythia8_AU2CT10_ggH400_ZZ4lep", kTRUE, 400, MC::ggF);
      m_sample[160181] = YieldSample("PowhegPythia8_AU2CT10_ggH420_ZZ4lep", kTRUE, 420, MC::ggF);
      m_sample[160182] = YieldSample("PowhegPythia8_AU2CT10_ggH440_ZZ4lep", kTRUE, 440, MC::ggF);
      m_sample[160183] = YieldSample("PowhegPythia8_AU2CT10_ggH460_ZZ4lep", kTRUE, 460, MC::ggF);
      m_sample[160184] = YieldSample("PowhegPythia8_AU2CT10_ggH480_ZZ4lep", kTRUE, 480, MC::ggF);
      m_sample[160185] = YieldSample("PowhegPythia8_AU2CT10_ggH500_ZZ4lep", kTRUE, 500, MC::ggF);
      m_sample[160186] = YieldSample("PowhegPythia8_AU2CT10_ggH520_ZZ4lep", kTRUE, 520, MC::ggF);
      m_sample[160187] = YieldSample("PowhegPythia8_AU2CT10_ggH540_ZZ4lep", kTRUE, 540, MC::ggF);
      m_sample[160188] = YieldSample("PowhegPythia8_AU2CT10_ggH560_ZZ4lep", kTRUE, 560, MC::ggF);
      m_sample[160189] = YieldSample("PowhegPythia8_AU2CT10_ggH580_ZZ4lep", kTRUE, 580, MC::ggF);
      m_sample[160190] = YieldSample("PowhegPythia8_AU2CT10_ggH600_ZZ4lep", kTRUE, 600, MC::ggF);
      m_sample[160202] = YieldSample("PowhegPythia8_AU2CT10_VBFH110_ZZ4lep", kTRUE, 110, MC::VBF);
      m_sample[160203] = YieldSample("PowhegPythia8_AU2CT10_VBFH115_ZZ4lep", kTRUE, 115, MC::VBF);
      m_sample[160204] = YieldSample("PowhegPythia8_AU2CT10_VBFH120_ZZ4lep", kTRUE, 120, MC::VBF);
      m_sample[160205] = YieldSample("PowhegPythia8_AU2CT10_VBFH125_ZZ4lep", kTRUE, 125, MC::VBF);
      m_sample[160206] = YieldSample("PowhegPythia8_AU2CT10_VBFH130_ZZ4lep", kTRUE, 130, MC::VBF);
      m_sample[160207] = YieldSample("PowhegPythia8_AU2CT10_VBFH135_ZZ4lep", kTRUE, 135, MC::VBF);
      m_sample[160208] = YieldSample("PowhegPythia8_AU2CT10_VBFH140_ZZ4lep", kTRUE, 140, MC::VBF);
      m_sample[160209] = YieldSample("PowhegPythia8_AU2CT10_VBFH145_ZZ4lep", kTRUE, 145, MC::VBF);
      m_sample[160210] = YieldSample("PowhegPythia8_AU2CT10_VBFH150_ZZ4lep", kTRUE, 150, MC::VBF);
      m_sample[160211] = YieldSample("PowhegPythia8_AU2CT10_VBFH155_ZZ4lep", kTRUE, 155, MC::VBF);
      m_sample[160212] = YieldSample("PowhegPythia8_AU2CT10_VBFH160_ZZ4lep", kTRUE, 160, MC::VBF);
      m_sample[160213] = YieldSample("PowhegPythia8_AU2CT10_VBFH165_ZZ4lep", kTRUE, 165, MC::VBF);
      m_sample[160214] = YieldSample("PowhegPythia8_AU2CT10_VBFH170_ZZ4lep", kTRUE, 170, MC::VBF);
      m_sample[160215] = YieldSample("PowhegPythia8_AU2CT10_VBFH175_ZZ4lep", kTRUE, 175, MC::VBF);
      m_sample[160216] = YieldSample("PowhegPythia8_AU2CT10_VBFH180_ZZ4lep", kTRUE, 180, MC::VBF);
      m_sample[160217] = YieldSample("PowhegPythia8_AU2CT10_VBFH185_ZZ4lep", kTRUE, 185, MC::VBF);
      m_sample[160218] = YieldSample("PowhegPythia8_AU2CT10_VBFH190_ZZ4lep", kTRUE, 190, MC::VBF);
      m_sample[160219] = YieldSample("PowhegPythia8_AU2CT10_VBFH195_ZZ4lep", kTRUE, 195, MC::VBF);
      m_sample[160220] = YieldSample("PowhegPythia8_AU2CT10_VBFH200_ZZ4lep", kTRUE, 200, MC::VBF);
      m_sample[160221] = YieldSample("PowhegPythia8_AU2CT10_VBFH220_ZZ4lep", kTRUE, 220, MC::VBF);
      m_sample[160222] = YieldSample("PowhegPythia8_AU2CT10_VBFH240_ZZ4lep", kTRUE, 240, MC::VBF);
      m_sample[160223] = YieldSample("PowhegPythia8_AU2CT10_VBFH260_ZZ4lep", kTRUE, 260, MC::VBF);
      m_sample[160224] = YieldSample("PowhegPythia8_AU2CT10_VBFH280_ZZ4lep", kTRUE, 280, MC::VBF);
      m_sample[160225] = YieldSample("PowhegPythia8_AU2CT10_VBFH300_ZZ4lep", kTRUE, 300, MC::VBF);
      m_sample[160226] = YieldSample("PowhegPythia8_AU2CT10_VBFH320_ZZ4lep", kTRUE, 320, MC::VBF);
      m_sample[160227] = YieldSample("PowhegPythia8_AU2CT10_VBFH340_ZZ4lep", kTRUE, 340, MC::VBF);
      m_sample[160228] = YieldSample("PowhegPythia8_AU2CT10_VBFH360_ZZ4lep", kTRUE, 360, MC::VBF);
      m_sample[160229] = YieldSample("PowhegPythia8_AU2CT10_VBFH380_ZZ4lep", kTRUE, 380, MC::VBF);
      m_sample[160230] = YieldSample("PowhegPythia8_AU2CT10_VBFH400_ZZ4lep", kTRUE, 400, MC::VBF);
      m_sample[160231] = YieldSample("PowhegPythia8_AU2CT10_VBFH420_ZZ4lep", kTRUE, 420, MC::VBF);
      m_sample[160232] = YieldSample("PowhegPythia8_AU2CT10_VBFH440_ZZ4lep", kTRUE, 440, MC::VBF);
      m_sample[160233] = YieldSample("PowhegPythia8_AU2CT10_VBFH460_ZZ4lep", kTRUE, 460, MC::VBF);
      m_sample[160234] = YieldSample("PowhegPythia8_AU2CT10_VBFH480_ZZ4lep", kTRUE, 480, MC::VBF);
      m_sample[160235] = YieldSample("PowhegPythia8_AU2CT10_VBFH500_ZZ4lep", kTRUE, 500, MC::VBF);
      m_sample[160236] = YieldSample("PowhegPythia8_AU2CT10_VBFH520_ZZ4lep", kTRUE, 520, MC::VBF);
      m_sample[160237] = YieldSample("PowhegPythia8_AU2CT10_VBFH540_ZZ4lep", kTRUE, 540, MC::VBF);
      m_sample[160238] = YieldSample("PowhegPythia8_AU2CT10_VBFH560_ZZ4lep", kTRUE, 560, MC::VBF);
      m_sample[160239] = YieldSample("PowhegPythia8_AU2CT10_VBFH580_ZZ4lep", kTRUE, 580, MC::VBF);
      m_sample[160240] = YieldSample("PowhegPythia8_AU2CT10_VBFH600_ZZ4lep", kTRUE, 600, MC::VBF);
      m_sample[160250] = YieldSample("Pythia8_AU2CTEQ6L1_WH100_ZZ4lep", kTRUE, 100, MC::WH);
      m_sample[160251] = YieldSample("Pythia8_AU2CTEQ6L1_WH105_ZZ4lep", kTRUE, 105, MC::WH);
      m_sample[160252] = YieldSample("Pythia8_AU2CTEQ6L1_WH110_ZZ4lep", kTRUE, 110, MC::WH);
      m_sample[160253] = YieldSample("Pythia8_AU2CTEQ6L1_WH115_ZZ4lep", kTRUE, 115, MC::WH);
      m_sample[160254] = YieldSample("Pythia8_AU2CTEQ6L1_WH120_ZZ4lep", kTRUE, 120, MC::WH);
      m_sample[160255] = YieldSample("Pythia8_AU2CTEQ6L1_WH125_ZZ4lep", kTRUE, 125, MC::WH);
      m_sample[160256] = YieldSample("Pythia8_AU2CTEQ6L1_WH130_ZZ4lep", kTRUE, 130, MC::WH);
      m_sample[160257] = YieldSample("Pythia8_AU2CTEQ6L1_WH135_ZZ4lep", kTRUE, 135, MC::WH);
      m_sample[160258] = YieldSample("Pythia8_AU2CTEQ6L1_WH140_ZZ4lep", kTRUE, 140, MC::WH);
      m_sample[160259] = YieldSample("Pythia8_AU2CTEQ6L1_WH145_ZZ4lep", kTRUE, 145, MC::WH);
      m_sample[160260] = YieldSample("Pythia8_AU2CTEQ6L1_WH150_ZZ4lep", kTRUE, 150, MC::WH);
      m_sample[160261] = YieldSample("Pythia8_AU2CTEQ6L1_WH155_ZZ4lep", kTRUE, 155, MC::WH);
      m_sample[160262] = YieldSample("Pythia8_AU2CTEQ6L1_WH160_ZZ4lep", kTRUE, 160, MC::WH);
      m_sample[160263] = YieldSample("Pythia8_AU2CTEQ6L1_WH165_ZZ4lep", kTRUE, 165, MC::WH);
      m_sample[160264] = YieldSample("Pythia8_AU2CTEQ6L1_WH170_ZZ4lep", kTRUE, 170, MC::WH);
      m_sample[160265] = YieldSample("Pythia8_AU2CTEQ6L1_WH175_ZZ4lep", kTRUE, 175, MC::WH);
      m_sample[160266] = YieldSample("Pythia8_AU2CTEQ6L1_WH180_ZZ4lep", kTRUE, 180, MC::WH);
      m_sample[160267] = YieldSample("Pythia8_AU2CTEQ6L1_WH185_ZZ4lep", kTRUE, 185, MC::WH);
      m_sample[160268] = YieldSample("Pythia8_AU2CTEQ6L1_WH190_ZZ4lep", kTRUE, 190, MC::WH);
      m_sample[160269] = YieldSample("Pythia8_AU2CTEQ6L1_WH195_ZZ4lep", kTRUE, 195, MC::WH);
      m_sample[160270] = YieldSample("Pythia8_AU2CTEQ6L1_WH200_ZZ4lep", kTRUE, 200, MC::WH);
      m_sample[160271] = YieldSample("Pythia8_AU2CTEQ6L1_WH220_ZZ4lep", kTRUE, 220, MC::WH);
      m_sample[160272] = YieldSample("Pythia8_AU2CTEQ6L1_WH240_ZZ4lep", kTRUE, 240, MC::WH);
      m_sample[160273] = YieldSample("Pythia8_AU2CTEQ6L1_WH260_ZZ4lep", kTRUE, 260, MC::WH);
      m_sample[160274] = YieldSample("Pythia8_AU2CTEQ6L1_WH280_ZZ4lep", kTRUE, 280, MC::WH);
      m_sample[160275] = YieldSample("Pythia8_AU2CTEQ6L1_WH300_ZZ4lep", kTRUE, 300, MC::WH);
      m_sample[160276] = YieldSample("Pythia8_AU2CTEQ6L1_WH320_ZZ4lep", kTRUE, 320, MC::WH);
      m_sample[160277] = YieldSample("Pythia8_AU2CTEQ6L1_WH340_ZZ4lep", kTRUE, 340, MC::WH);
      m_sample[160278] = YieldSample("Pythia8_AU2CTEQ6L1_WH360_ZZ4lep", kTRUE, 360, MC::WH);
      m_sample[160279] = YieldSample("Pythia8_AU2CTEQ6L1_WH380_ZZ4lep", kTRUE, 380, MC::WH);
      m_sample[160280] = YieldSample("Pythia8_AU2CTEQ6L1_WH400_ZZ4lep", kTRUE, 400, MC::WH);
      m_sample[160300] = YieldSample("Pythia8_AU2CTEQ6L1_ZH100_ZZ4lep", kTRUE, 100, MC::ZH);
      m_sample[160301] = YieldSample("Pythia8_AU2CTEQ6L1_ZH105_ZZ4lep", kTRUE, 105, MC::ZH);
      m_sample[160302] = YieldSample("Pythia8_AU2CTEQ6L1_ZH110_ZZ4lep", kTRUE, 110, MC::ZH);
      m_sample[160303] = YieldSample("Pythia8_AU2CTEQ6L1_ZH115_ZZ4lep", kTRUE, 115, MC::ZH);
      m_sample[160304] = YieldSample("Pythia8_AU2CTEQ6L1_ZH120_ZZ4lep", kTRUE, 120, MC::ZH);
      m_sample[160305] = YieldSample("Pythia8_AU2CTEQ6L1_ZH125_ZZ4lep", kTRUE, 125, MC::ZH);
      m_sample[160306] = YieldSample("Pythia8_AU2CTEQ6L1_ZH130_ZZ4lep", kTRUE, 130, MC::ZH);
      m_sample[160307] = YieldSample("Pythia8_AU2CTEQ6L1_ZH135_ZZ4lep", kTRUE, 135, MC::ZH);
      m_sample[160308] = YieldSample("Pythia8_AU2CTEQ6L1_ZH140_ZZ4lep", kTRUE, 140, MC::ZH);
      m_sample[160309] = YieldSample("Pythia8_AU2CTEQ6L1_ZH145_ZZ4lep", kTRUE, 145, MC::ZH);
      m_sample[160310] = YieldSample("Pythia8_AU2CTEQ6L1_ZH150_ZZ4lep", kTRUE, 150, MC::ZH);
      m_sample[160311] = YieldSample("Pythia8_AU2CTEQ6L1_ZH155_ZZ4lep", kTRUE, 155, MC::ZH);
      m_sample[160312] = YieldSample("Pythia8_AU2CTEQ6L1_ZH160_ZZ4lep", kTRUE, 160, MC::ZH);
      m_sample[160313] = YieldSample("Pythia8_AU2CTEQ6L1_ZH165_ZZ4lep", kTRUE, 165, MC::ZH);
      m_sample[160314] = YieldSample("Pythia8_AU2CTEQ6L1_ZH170_ZZ4lep", kTRUE, 170, MC::ZH);
      m_sample[160315] = YieldSample("Pythia8_AU2CTEQ6L1_ZH175_ZZ4lep", kTRUE, 175, MC::ZH);
      m_sample[160316] = YieldSample("Pythia8_AU2CTEQ6L1_ZH180_ZZ4lep", kTRUE, 180, MC::ZH);
      m_sample[160317] = YieldSample("Pythia8_AU2CTEQ6L1_ZH185_ZZ4lep", kTRUE, 185, MC::ZH);
      m_sample[160318] = YieldSample("Pythia8_AU2CTEQ6L1_ZH190_ZZ4lep", kTRUE, 190, MC::ZH);
      m_sample[160319] = YieldSample("Pythia8_AU2CTEQ6L1_ZH195_ZZ4lep", kTRUE, 195, MC::ZH);
      m_sample[160320] = YieldSample("Pythia8_AU2CTEQ6L1_ZH200_ZZ4lep", kTRUE, 200, MC::ZH);
      m_sample[160321] = YieldSample("Pythia8_AU2CTEQ6L1_ZH220_ZZ4lep", kTRUE, 220, MC::ZH);
      m_sample[160322] = YieldSample("Pythia8_AU2CTEQ6L1_ZH240_ZZ4lep", kTRUE, 240, MC::ZH);
      m_sample[160323] = YieldSample("Pythia8_AU2CTEQ6L1_ZH260_ZZ4lep", kTRUE, 260, MC::ZH);
      m_sample[160324] = YieldSample("Pythia8_AU2CTEQ6L1_ZH280_ZZ4lep", kTRUE, 280, MC::ZH);
      m_sample[160325] = YieldSample("Pythia8_AU2CTEQ6L1_ZH300_ZZ4lep", kTRUE, 300, MC::ZH);
      m_sample[160326] = YieldSample("Pythia8_AU2CTEQ6L1_ZH320_ZZ4lep", kTRUE, 320, MC::ZH);
      m_sample[160327] = YieldSample("Pythia8_AU2CTEQ6L1_ZH340_ZZ4lep", kTRUE, 340, MC::ZH);
      m_sample[160328] = YieldSample("Pythia8_AU2CTEQ6L1_ZH360_ZZ4lep", kTRUE, 360, MC::ZH);
      m_sample[160329] = YieldSample("Pythia8_AU2CTEQ6L1_ZH380_ZZ4lep", kTRUE, 380, MC::ZH);
      m_sample[160330] = YieldSample("Pythia8_AU2CTEQ6L1_ZH400_ZZ4lep", kTRUE, 400, MC::ZH);
      m_sample[146980] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZeeNp0", kFALSE, -1, MC::none);
      m_sample[146981] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZeeNp1", kFALSE, -1, MC::none);
      m_sample[146982] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZbbeeNp2", kFALSE, -1, MC::none);
      m_sample[146985] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZmumuNp0", kFALSE, -1, MC::none);
      m_sample[146986] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZmumuNp1", kFALSE, -1, MC::none);
      m_sample[146987] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_4lFilter_ZbbmumuNp2", kFALSE, -1, MC::none);
      m_sample[146990] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZeeNp0", kFALSE, -1, MC::none);
      m_sample[146991] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZeeNp1", kFALSE, -1, MC::none);
      m_sample[146992] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZbbeeNp2", kFALSE, -1, MC::none);
      m_sample[146995] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZmumuNp0", kFALSE, -1, MC::none);
      m_sample[146996] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZmumuNp1", kFALSE, -1, MC::none);
      m_sample[146997] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_3lFilter_4lVeto_ZbbmumuNp2", kFALSE, -1, MC::none);
      m_sample[107650] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp0", kFALSE, -1, MC::Z);
      m_sample[107651] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp1", kFALSE, -1, MC::Z);
      m_sample[107652] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp2", kFALSE, -1, MC::Z);
      m_sample[107653] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp3", kFALSE, -1, MC::Z);
      m_sample[107654] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp4", kFALSE, -1, MC::Z);
      m_sample[107655] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZeeNp5", kFALSE, -1, MC::Z);
      m_sample[107660] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp0", kFALSE, -1, MC::Z);
      m_sample[107661] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp1", kFALSE, -1, MC::Z);
      m_sample[107662] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp2", kFALSE, -1, MC::Z);
      m_sample[107663] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp3", kFALSE, -1, MC::Z);
      m_sample[107664] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp4", kFALSE, -1, MC::Z);
      m_sample[107665] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp5", kFALSE, -1, MC::Z);
      m_sample[107670] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp0", kFALSE, -1, MC::Z);
      m_sample[107671] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp1", kFALSE, -1, MC::Z);
      m_sample[107672] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp2", kFALSE, -1, MC::Z);
      m_sample[107673] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp3", kFALSE, -1, MC::Z);
      m_sample[107674] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp4", kFALSE, -1, MC::Z);
      m_sample[107675] = YieldSample("AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp5", kFALSE, -1, MC::Z);
      m_sample[146830] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp0Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146831] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp1Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146832] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp2Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146833] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp3Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146834] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp4Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146840] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp0Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146841] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp1Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146842] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp2Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146843] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp3Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146844] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp4Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146850] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp0Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146851] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp1Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146852] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp2Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[146853] = YieldSample("AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp3Excl_Mll10to60", kFALSE, -1, MC::Z);
      m_sample[147806] = YieldSample("PowhegPythia8_AU2CT10_Zee", kFALSE, -1, MC::none);
      m_sample[147806] = YieldSample("PowhegPythia8_AU2CT10_Zee", kFALSE, -1, MC::none);
      m_sample[147807] = YieldSample("PowhegPythia8_AU2CT10_Zmumu", kFALSE, -1, MC::none);
      m_sample[126937] = YieldSample("PowhegPythia8_AU2CT10_ZZ_4e_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[126938] = YieldSample("PowhegPythia8_AU2CT10_ZZ_2e2mu_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[126939] = YieldSample("PowhegPythia8_AU2CT10_ZZ_2e2tau_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[126940] = YieldSample("PowhegPythia8_AU2CT10_ZZ_4mu_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[126941] = YieldSample("PowhegPythia8_AU2CT10_ZZ_2mu2tau_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[126942] = YieldSample("PowhegPythia8_AU2CT10_ZZ_4tau_mll4_2pt5", kFALSE, -1, MC::ZZ);
      m_sample[116600] = YieldSample("gg2ZZJimmy_AUET2CT10_ZZ4lep", kFALSE, -1, MC::none);
      m_sample[116601] = YieldSample("gg2ZZJimmy_AUET2CT10_ZZ4e", kFALSE, -1, MC::gg2ZZ);
      m_sample[116602] = YieldSample("gg2ZZJimmy_AUET2CT10_ZZ4mu", kFALSE, -1, MC::gg2ZZ);
      m_sample[116603] = YieldSample("gg2ZZJimmy_AUET2CT10_ZZ2e2mu", kFALSE, -1, MC::gg2ZZ);
      m_sample[126894] = YieldSample("Sherpa_CT10_llll_ZZ", kFALSE, -1, MC::none);
      m_sample[105200] = YieldSample("McAtNloJimmy_CT10_ttbar_LeptonFilter", kFALSE, -1, MC::tt);
      m_sample[110001] = YieldSample("McAtNloJimmy_CT10_ttbar_dilepton", kFALSE, -1, MC::none);
      m_sample[146369] = YieldSample("McAtNloJimmy_CT10_ttbar_4LepMass_Mll40GeV12GeV", kFALSE, -1, MC::tt);

      std::cout << "INFO : filled map of mc12a samples" << std::endl;
      return kTRUE;
   }

   return kFALSE;
}

void YieldEstimator::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   fillMCMap("mc12a");
   m_COM_energy = CrossSections::EightTeV;
   lumi[0] = 5.831; // 4mu
   lumi[1] = 5.831; // 2mu2e
   lumi[2] = 5.831; // 2e2mu
   lumi[3] = 5.858; // 4e

   m_output_raw = new TFile("output_raw.root", "RECREATE");
   m_output_with_constraint = new TFile("output_constraint.root", "RECREATE");
   m_output_without_constraint = new TFile("output_no_constraint.root", "RECREATE");

   m_histos["constr_1300"] = FinalHistograms(0.5);
   m_histos["constr_1300"].setDirectory(m_output_with_constraint);

}

void YieldEstimator::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t YieldEstimator::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either YieldEstimator::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   if (!m_interesting) return kTRUE;

   fChain->GetTree()->GetEntry(entry);

   if (run != m_current_run) {
      Error("Process", "read event with run number %d which is not the same as %d got from the first entry: everything will screw up!", run, m_current_run);
   }

   Double_t poids(1.);

   poids = pu_weight * trigSF_weight * Z1_lepplus_weight * Z1_lepminus_weight * Z2_lepplus_weight * Z2_lepminus_weight * top_weight * powhegbug_weight * vxz_weight;
   poids *= lumi[type] * m_sample[run].xsec;// / generated[run];

   if (selected == 1) {
      m_sample[run].h_H_m[type]->Fill(H_m / 1000., poids);
      m_sample[run].h_H_m_constrained[type]->Fill(H_m_constrained / 1000., poids);
      m_sample[run].h_Z1_m[type]->Fill(Z1_m / 1000., poids);
      m_sample[run].h_Z2_m[type]->Fill(Z2_m / 1000., poids);
   }


   return kTRUE;
}

void YieldEstimator::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   // complete the rescale to luminosity
   for (std::map<UInt_t, YieldSample>::iterator sample = m_sample.begin(); sample != m_sample.end(); ++sample) {
      const UInt_t thisRun = sample->first;
      const YieldSample &thisSample = sample->second;

      if (thisSample.generated > 0)
         thisSample.print();

      const Float_t one_over_generated = 1. / thisSample.generated;

      for (UInt_t final_state = 0; final_state < thisSample.h_H_m.size(); final_state++) {
         m_sample[thisRun].h_H_m[final_state]->Scale(one_over_generated);
         m_sample[thisRun].h_H_m_constrained[final_state]->Scale(one_over_generated);
         m_sample[thisRun].h_Z1_m[final_state]->Scale(one_over_generated);
         m_sample[thisRun].h_Z2_m[final_state]->Scale(one_over_generated);
         // combine into the various histos (todo: use a map)
         if (m_sample[thisRun].type == MC::ggF) {
            if (m_sample[thisRun].mass == 130)
               m_histos["constr_1300"].higgs[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::VBF) {
            if (m_sample[thisRun].mass == 130)
               m_histos["constr_1300"].higgsVBF[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::WH) {
            if (m_sample[thisRun].mass == 130)
               m_histos["constr_1300"].higgsWH[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::ZH) {
            if (m_sample[thisRun].mass == 130)
               m_histos["constr_1300"].higgsZH[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::Z) {
            m_histos["constr_1300"].histoZ[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::Zbb) {
            m_histos["constr_1300"].histoZbb[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::tt) {
            m_histos["constr_1300"].histott[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::ZZ) {
            m_histos["constr_1300"].histoZZ[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         } else if (m_sample[thisRun].type == MC::gg2ZZ) {
            m_histos["constr_1300"].histogg2ZZ[final_state]->Add(m_sample[thisRun].h_H_m[final_state]);
         }
      }
   }

}

void YieldEstimator::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   m_output_raw->Write();
   m_output_raw->Close();
   m_output_with_constraint->Write();
   m_output_with_constraint->Close();
   m_output_without_constraint->Write();
   m_output_without_constraint->Close();
}
