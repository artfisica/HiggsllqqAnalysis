from ROOT    import *
from utils   import *
from samples import *

############################################################# 2012 SAMPLES ###########################################################

Top               = extendedsample()
Top.name          = 'Top'
Top.color         = kOrange-3
Top.isMC          = True
Top.isBKG         = True
Top.samples.append(McAtNloJimmy_CT10_ttbar_LeptonFilter)
Top.samples.append(McAtNloJimmy_AUET2CT10_SingleTopWtChanIncl)
Top.samples.append(AcerMCPythia_AUET2BCTEQ6L1_singletop_tchan_e)
Top.samples.append(AcerMCPythia_AUET2BCTEQ6L1_singletop_tchan_mu)
Top.samples.append(AcerMCPythia_AUET2BCTEQ6L1_singletop_tchan_tau)


Diboson           = extendedsample()
Diboson.name      = 'Diboson'
Diboson.color     = kGreen
Diboson.isMC      = True
Diboson.isBKG     = True
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmenuenu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmenumunu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmenutaunu) 
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmmunuenu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmmunumunu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmmunutaunu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmtaunuenu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmtaunumunu)
## Diboson.samples.append(McAtNloJimmy_AUET2CT10_WpWmtaunutaunu) 
Diboson.samples.append(Herwig_AUET2CTEQ6L1_WW)
Diboson.samples.append(Herwig_AUET2CTEQ6L1_ZZ)
Diboson.samples.append(Herwig_AUET2CTEQ6L1_WZ)


Zlf_Alp            = extendedsample()
Zlf_Alp.name       = 'Zlf_Alp'
Zlf_Alp.color      = kMagenta
Zlf_Alp.isMC       = True
Zlf_Alp.isBKG      = True
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp4)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp5)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp4)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeNp5)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp4)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp5)
#### DY
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp0Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp1Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp2Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp3Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp4Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp5Incl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp0Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp1Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp2Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp3Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp4Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp5Incl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp0Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp1Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp2Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp3Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp4Excl_Mll10to60)
Zlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp5Incl_Mll10to60)
#### Z HF
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp3)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp0)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp1)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp2)
Zlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp3)






DYlf_Alp           = extendedsample()
DYlf_Alp.name      = 'DYlf_Alp'
DYlf_Alp.color     = kMagenta+2
DYlf_Alp.isMC      = True
DYlf_Alp.isBKG     = True
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp0Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp1Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp2Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp3Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp4Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp5Incl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp0Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp1Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp2Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp3Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp4Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp5Incl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp0Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp1Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp2Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp3Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp4Excl_Mll10to60)
DYlf_Alp.samples.append(AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp5Incl_Mll10to60)






ZHf_Alp            = extendedsample()
ZHf_Alp.name       = 'ZHf_Alp'
ZHf_Alp.color      = kMagenta-2
ZHf_Alp.isMC       = True
ZHf_Alp.isBKG      = True
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautauccNp3)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeeccNp3)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumuccNp3)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZeebbNp3)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZmumubbNp3)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp0)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp1)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp2)
ZHf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_ZtautaubbNp3)




Wlf_Alp           = extendedsample()
Wlf_Alp.name      = 'Wlf_Alp'
Wlf_Alp.color     = kMagenta-2
Wlf_Alp.isMC      = True
Wlf_Alp.isBKG     = True
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WbbNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WbbNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WbbNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WbbNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp4)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WenuNp5)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp4)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WmunuNp5)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp4)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp5)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WccNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WccNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WccNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WccNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WcNp0)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WcNp1)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WcNp2)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WcNp3)
Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WcNp4)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp0)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp1)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp2)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp3)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp4)
## Wlf_Alp.samples.append(AlpgenJimmy_AUET2CTEQ6L1_WgammaNp5)









Signal_ggH        = extendedsample()
Signal_ggH.name   = 'Signal_ggH'
Signal_ggH.color  = kRed
Signal_ggH.isMC   = True
Signal_ggH.isBKG  = False
## Go into plotterino.py to choose the Signal sample to use


Signal_VBFH       = extendedsample()
Signal_VBFH.name  = 'Signal_VBFH'
Signal_VBFH.color = kRed
Signal_VBFH.isMC  = True
Signal_VBFH.isBKG = False
## Go into plotterino.py to choose the Signal sample to use


Data              = extendedsample()
Data.name         = 'Data'
Data.color        = kBlack
Data.isMC         = False
Data.isBKG        = False
Data.samples.append(physics_Egamma)
Data.samples.append(physics_Muons)


QCD               = extendedsample()
QCD.name          = 'QCD'
QCD.color         = kOrange+4
QCD.isMC          = False
QCD.isBKG         = True
QCD.samples.append(physics_Egamma)
QCD.samples.append(physics_Muons)


TotalBkg          = extendedsample()
TotalBkg.name     = 'totbkg'
TotalBkg.color    = kBlack
TotalBkg.isMC     = True
TotalBkg.isBKG    = True
