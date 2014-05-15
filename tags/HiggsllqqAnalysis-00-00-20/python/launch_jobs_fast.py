#!/bin/bash

data = [
  'user.vippolit.ZNtuple000208.run_177986_191933.physics_Muons.rootcore_r17.01/',
  'user.vippolit.ZNtuple000208.run_177986_191933.physics_Egamma.rootcore_r17.01/',
]

AlpgenJimmyW_Samples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.107693.AlpgenJimmyWmunuNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107694.AlpgenJimmyWmunuNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107690.AlpgenJimmyWmunuNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107692.AlpgenJimmyWmunuNp2_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107684.AlpgenJimmyWenuNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107701.AlpgenJimmyWtaunuNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107700.AlpgenJimmyWtaunuNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107705.AlpgenJimmyWtaunuNp5_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107691.AlpgenJimmyWmunuNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107683.AlpgenJimmyWenuNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107685.AlpgenJimmyWenuNp5_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107703.AlpgenJimmyWtaunuNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107682.AlpgenJimmyWenuNp2_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107681.AlpgenJimmyWenuNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107704.AlpgenJimmyWtaunuNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107680.AlpgenJimmyWenuNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107702.AlpgenJimmyWtaunuNp2_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107695.AlpgenJimmyWmunuNp5_pt20.merge.NTUP_HSG2..qqll.01/",
  ###
  "user.arturos.ZNtuple000207.mc11_7TeV.106043.PythiaWenu_no_filter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.106044.PythiaWmunu_no_filter.merge.NTUP_HSG2..qqll.01/",
]

AlpgenJimmyZ_Samples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.107652.AlpgenJimmyZeeNp2_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107655.AlpgenJimmyZeeNp5_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107654.AlpgenJimmyZeeNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107653.AlpgenJimmyZeeNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107650.AlpgenJimmyZeeNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107651.AlpgenJimmyZeeNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107661.AlpgenJimmyZmumuNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107662.AlpgenJimmyZmumuNp2_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107663.AlpgenJimmyZmumuNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107665.AlpgenJimmyZmumuNp5_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107660.AlpgenJimmyZmumuNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107664.AlpgenJimmyZmumuNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107664.AlpgenJimmyZmumuNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107670.AlpgenJimmyZtautauNp0_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107671.AlpgenJimmyZtautauNp1_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107673.AlpgenJimmyZtautauNp3_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107674.AlpgenJimmyZtautauNp4_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.107675.AlpgenJimmyZtautauNp5_pt20.merge.NTUP_HSG2..qqll.01/",
]

AlpgenJimmyZbb_Samples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.109308.AlpgenJimmyZmumubbNp3_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109307.AlpgenJimmyZmumubbNp2_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109305.AlpgenJimmyZmumubbNp0_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109306.AlpgenJimmyZmumubbNp1_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109310.AlpgenJimmyZtautaubbNp0_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109311.AlpgenJimmyZtautaubbNp1_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109312.AlpgenJimmyZtautaubbNp2_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109313.AlpgenJimmyZtautaubbNp3_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109300.AlpgenJimmyZeebbNp0_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109302.AlpgenJimmyZeebbNp2_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109301.AlpgenJimmyZeebbNp1_nofilter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109303.AlpgenJimmyZeebbNp3_nofilter.merge.NTUP_HSG2..qqll.01/",
]

DrellYan_Samples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.116251.AlpgenJimmyZeeNp1_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116270.AlpgenJimmyZtautauNp0_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116262.AlpgenJimmyZmumuNp2_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116260.AlpgenJimmyZmumuNp0_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116255.AlpgenJimmyZeeNp5_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116275.AlpgenJimmyZtautauNp5_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116264.AlpgenJimmyZmumuNp4_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116252.AlpgenJimmyZeeNp2_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116274.AlpgenJimmyZtautauNp4_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116261.AlpgenJimmyZmumuNp1_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116272.AlpgenJimmyZtautauNp2_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116273.AlpgenJimmyZtautauNp3_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116254.AlpgenJimmyZeeNp4_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116250.AlpgenJimmyZeeNp0_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116265.AlpgenJimmyZmumuNp5_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116271.AlpgenJimmyZtautauNp1_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116253.AlpgenJimmyZeeNp3_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116263.AlpgenJimmyZmumuNp3_Mll10to40_pt20.merge.NTUP_HSG2..qqll.01/",
  ###
  "user.arturos.ZNtuple000207.mc11_7TeV.108322.PythiaDrellYanLowM_ee3.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108321.PythiaDrellYanLowM_mu3.merge.NTUP_HSG2..qqll.01/",
  ###
  "user.arturos.ZNtuple000207.mc11_7TeV.128813.SherpaZZllll.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116600.gg2ZZ_JIMMY_ZZ4lep.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126406.PowHegBoxZZeenn_Pythia.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126407.PowHegBoxZZmmnn_Pythia.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126408.PowHegBoxZZttnn_Pythia.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126401.PowHegBoxZZmmmm_Pythia_mll025_m4l40.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126399.PowHegBoxZZeemm_Pythia_mll025_m4l40.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126400.PowHegBoxZZeeee_Pythia_mll025_m4l40.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126417.AlpgenJimmyZeeccNp3_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126416.AlpgenJimmyZeeccNp2_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126414.AlpgenJimmyZeeccNp0_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126415.AlpgenJimmyZeeccNp1_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126419.AlpgenJimmyZmumuccNp1_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126421.AlpgenJimmyZmumuccNp3_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126418.AlpgenJimmyZmumuccNp0_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.126420.AlpgenJimmyZmumuccNp2_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128135.AlpgenJimmyLowMassDYmumubbNp0_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128136.AlpgenJimmyLowMassDYmumubbNp1_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128138.AlpgenJimmyLowMassDYmumubbNp3_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128137.AlpgenJimmyLowMassDYmumubbNp2_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128130.AlpgenJimmyLowMassDYeebbNp0_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128131.AlpgenJimmyLowMassDYeebbNp1_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128132.AlpgenJimmyLowMassDYeebbNp2_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128133.AlpgenJimmyLowMassDYeebbNp3_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128141.AlpgenJimmyLowMassDYtautaubbNp1_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128142.AlpgenJimmyLowMassDYtautaubbNp2_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128140.AlpgenJimmyLowMassDYtautaubbNp0_nofilter.merge.NTUP_HSG2.e.qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.128143.AlpgenJimmyLowMassDYtautaubbNp3_nofilter.merge.NTUP_HSG2.e.qqll.01/",
]

PowHegPythia_ggH_SignalSamples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.116830.PowHegPythia_ggH130_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116828.PowHegPythia_ggH120_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116841.PowHegPythia_ggH185_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116826.PowHegPythia_ggH110_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116857.PowHegPythia_ggH440_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116840.PowHegPythia_ggH180_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116843.PowHegPythia_ggH195_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116847.PowHegPythia_ggH220_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116838.PowHegPythia_ggH170_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116856.PowHegPythia_ggH420_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116858.PowHegPythia_ggH460_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116842.PowHegPythia_ggH190_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116854.PowHegPythia_ggH360_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116859.PowHegPythia_ggH480_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116834.PowHegPythia_ggH150_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116855.PowHegPythia_ggH380_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116852.PowHegPythia_ggH320_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116848.PowHegPythia_ggH240_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116863.PowHegPythia_ggH560_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116851.PowHegPythia_ggH300_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116850.PowHegPythia_ggH280_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116618.PowHegPythia_ggH400_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116844.PowHegPythia_ggH200_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116846.PowHegPythia_ggH210_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116853.PowHegPythia_ggH340_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116865.PowHegPythia_ggH600_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116864.PowHegPythia_ggH580_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116849.PowHegPythia_ggH260_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116832.PowHegPythia_ggH140_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116860.PowHegPythia_ggH500_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116836.PowHegPythia_ggH160_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116862.PowHegPythia_ggH540_ZZllqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.116861.PowHegPythia_ggH520_ZZllqq.merge.NTUP_HSG2..qqll.01/",
]

McAtNlo_JIMMY_Diboson_Samples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.105926.McAtNlo_JIMMY_WpWm_munutaunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105928.McAtNlo_JIMMY_WpWm_taunuenu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105922.McAtNlo_JIMMY_WpWm_enumunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105929.McAtNlo_JIMMY_WpWm_taunumunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105921.McAtNlo_JIMMY_WpWm_enuenu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105923.McAtNlo_JIMMY_WpWm_enutaunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105927.McAtNlo_JIMMY_WpWm_taunutaunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105924.McAtNlo_JIMMY_WpWm_munumunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105925.McAtNlo_JIMMY_WpWm_munuenu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105942.McAtNlo_JIMMY_WpZ_qqll.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105970.McAtNlo_JIMMY_WmZ_lnuqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105941.McAtNlo_JIMMY_WpZ_lnull.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105940.McAtNlo_JIMMY_WpZ_lnuqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105971.McAtNlo_JIMMY_WmZ_lnull.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105972.McAtNlo_JIMMY_WmZ_qqll.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105930.McAtNlo_JIMMY_ZZ_llqq.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105932.McAtNlo_JIMMY_ZZ_llnunu.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.105931.McAtNlo_JIMMY_ZZ_llll.merge.NTUP_HSG2..qqll.01/",
  ###
  "user.arturos.ZNtuple000207.mc11_7TeV.109963.McAtNlo_JIMMY_H110_WpWm_taunutaunu.merge.NTUP_HSG2..qqll.01/",
"user.arturos.ZNtuple000207.mc11_7TeV.109967.McAtNlo_JIMMY_H130_WpWm_taunutaunu.merge.NTUP_HSG2..qqll.01/",
]

TopSamples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.105200.T1_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108346.st_Wt_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108342.st_tchan_taunu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108345.st_schan_taunu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108340.st_tchan_enu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108343.st_schan_enu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.108344.st_schan_munu_McAtNlo_Jimmy.merge.NTUP_HSG2..qqll.01/",
]

PythiaSamples = [
  "user.arturos.ZNtuple000207.mc11_7TeV.109835.PythiaH150zzqqll.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109833.PythiaH400zzqqll.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.109360.PythiaH600zzqqll.merge.NTUP_HSG2..qqll.01/",
  ###
  "user.arturos.ZNtuple000207.mc11_7TeV.106046.PythiaZee_no_filter.merge.NTUP_HSG2..qqll.01/",
  "user.arturos.ZNtuple000207.mc11_7TeV.106047.PythiaZmumu_no_filter.merge.NTUP_HSG2..qqll.01/",
]


lists_to_process = [
  PythiaSamples
]


### analysis options
isMC           = True
CaloMuons      = False       # General -noSmearing Option
muonFamily     = "STACO"     # STACO or MUID or CALO 
electronFamily = "GSF"
analyses = [
  "MU2",
  "E2",
  ##"MUE"
]

#### Kind of Production #####
TypeOfProduction = "Private"          ## Change "Private" to "Official" JUST when you have been asked to make an official Hqqll Ntuple production!!!
#############################         ## For 28th March 2012 ==> ZNtuple000207.mc11_7TeV

### grid options
username="arturos"                    ## Who is running the job <username>
DoLowMass = True                      ## Change it propertly

if(TypeOfProduction == "Private"):
  mycodeversion = "myHqqll"           ## your private Conversion
  version       = "xx"                ## your private Version
if(TypeOfProduction == "Official"):
  mycodeversion = "officialHqqll"     ## Official Conversion from 28th March 2012. Please dont change it.
  version       = "00"                ## Official Version from 28th March 2012.
if(isMC):
  suffix= "MC."   + version             
else:
  suffix= "Data." + version 
  
filename="output.root"
excludedSite=""                       ## Example: "ANALY_AGLT2"

useItCloud   = False                  ## To use the Italian Cloud Only!
SendToNaples = False                  ## To commit the output files directly to Local disk in Naples TR2
SendToRome   = True                   ## To commit the output files directly to Local disk in Rome   TR2

if(isMC):
  isPlain = True                      ## Dataset name appended as is
  group   = False                     ## Run a single job over whole dataset_list
else:
  isPlain = False
  group   = False                     ## Run a single job over whole dataset_list
  
avoidBigJobs = True                   ## 10 Gb limit
period_name  = "run_177986_191933"
printOnly    = True                   ## Do not execute the command, just print it


# prun usually skips files you need! check its output and put them here
files_to_preserve = [
  'HiggsAnalysis/packages/files/pileup/*root*',
  'HiggsAnalysis/packages/files/ggFHiggsPtWeight/*root*',
  '/JetResolution/share/JERProviderPlots.root',
]

useTarBall = True                     ## To speed to sending of hte Jobs, keep it 'True'.

#####################
#####################

import os
import re

primo = True
for channel in analyses:
  _suffix = channel + "." + suffix

  for dataset_list in lists_to_process:
  
    if (not isMC): # real data
      runnumber = re.compile('\w*\.\w*.\w*\.(\d+)') 
      stream    = re.compile('physics_(\w*)') 
    else: # MC
      runnumber = re.compile('\w*\.\w*.\w*.(\w+)') # using i.e. mc09_7TeV
      stream    = re.compile('\w*\.\w*.\w*\.\w+\.\d+\.(\w+)')
    
    option_string = ""
    if (isMC):
      option_string += " -isMC"
    option_string += " -type " + channel
    option_string += " -input input.txt -output " + filename
    option_string += " -muonFamily " + muonFamily
    option_string += " -electronFamily " + electronFamily
    if (CaloMuons):
      option_string += " -noSmearing"
    if(DoLowMass):
      option_string += " -dolowmass"
      
    prun_exec_string = "--exec \"echo %IN > input.txt; python HiggsqqllAnalysis/python/regolari.py; cat input.txt; HiggsqqllAnalysis/bin/runAnalysis" + option_string + "\""
    
    
    if (group):
      dataset_string = ""
      for this_dataset in dataset_list:
        dataset_string += this_dataset + ","
      dataset_string = dataset_string[:-1]
    
      this_runnumber = period_name
      this_stream = stream.search(dataset_list[0]).group(1)
    
      if (not isMC):
        this_stream = "physics_" + this_stream
    
      this_command = "prun " + prun_exec_string + " --outDS user."+username+"." + mycodeversion + "." + this_runnumber + "." + this_stream + "." + _suffix + " --outputs " + filename + " --useAthenaPackages --inDS " + dataset_string
    
      if (avoidBigJobs):
        this_command += " --nGBPerJob=10"
    
      ## add the excluded sites
      if (excludedSite != ""):
        this_command += " --excludedSite " + excludedSite
  
      # add the rootcore stuff
      this_command += " --useRootCore"

      # preserve files
      if (len(files_to_preserve) > 0):
        this_command += " --extFile="
        for interesting_file in files_to_preserve:
          this_command += interesting_file + ','
      
      # use it cloud?
      if (useItCloud):
        this_command += " --cloud=IT"
    
      print this_command
      if (not printOnly):
        os.system(this_command)
    
    else:
      for this_dataset in dataset_list:
        if (isPlain):
          this_stream = this_dataset
          if (this_stream[-1:] == '/'):
            this_stream = this_stream[:-1]
          this_stream = this_stream[+37:][:-25]
          this_command = "prun " + prun_exec_string + " --outDS user."+username+"." + mycodeversion + "." + this_stream + "." + _suffix + " --outputs " + filename + " --useAthenaPackages --inDS " + this_dataset
        else:
          if (not isMC):
            this_runnumber = "run_" + runnumber.search(this_dataset).group(1)
          else:
            this_runnumber = runnumber.search(this_dataset).group(1)
          this_stream = stream.search(this_dataset).group(1)
    
          if (not isMC):
            this_stream = "physics_" + this_stream
          
          this_command = "prun " + prun_exec_string + " --outDS user."+username+"." + mycodeversion + "." + this_runnumber + "." + this_stream + "." + _suffix + " --outputs " + filename + " --useAthenaPackages --inDS " + this_dataset
    
        if (avoidBigJobs):
          this_command += " --nGBPerJob=10"
    
         ## add the excluded sites
        if (excludedSite != ""):
          this_command += " --excludedSite " + excludedSite
  
        # add the rootcore stuff
        this_command += " --useRootCore"

        if(SendToNaples):
          this_command += " --destSE INFN-NAPOLI-ATLAS_LOCALGROUPDISK"
        if(SendToRome):
          this_command += " --destSE INFN-ROMA1_LOCALGROUPDISK"

        # preserve files
        if (len(files_to_preserve) > 0):
          this_command += " --extFile="
          for interesting_file in files_to_preserve:
            this_command += interesting_file + ','

        # use it cloud?
        if (useItCloud):
          this_command += " --cloud=IT"
          
        primo_stream =""  
        if (useTarBall):
          if(primo):
            primo_stream += this_command
            primo_stream += " --outTarBall=tarball.tar --noSubmit"
          this_command += " --inTarBall=tarball.tar"

        if(useTarBall and primo):
          print primo_stream
          if (not printOnly):
            os.system(primo_stream)
        primo = False  
        print this_command
        
        if (not printOnly):
          os.system(this_command)

## Deleting the .tar file!
if (not printOnly):
  os.system("mv tarball.tar /tmp/"+username+"/.")
