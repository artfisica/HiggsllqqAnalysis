## October 7th 2013
dq2-ls -L INFN-ROMA1_LOCALGROUPDISK -fp  user.arturos.mc12_8TeV.117293.AlpgenJimmy_AUET2CTEQ6L1_WcNp0.merge.NTUP_HSG2.e1601_s1499_s1504_r3658_r3549_p1344.13.1/|grep srm|awk '{print "lcg-la "$1}'|awk '{system($0 "|sed 's#lfn:/grid#root://t2-dpm-01.na.infn.it:1094/#g'")}'

## November 28th 2013
 dq2-ls -L INFN-NAPOLI-ATLAS_PHYS-HIGGS -fp  mc12_8TeV.105200.McAtNloJimmy_CT10_ttbar_LeptonFilter.merge.NTUP_HSG2.e1513_s1499_s1504_r3658_r3549_p1344/|grep srm|awk '{print "lcg-la "$1}'|awk '{system($0 "|sed 's#lfn:/grid#root://t2-dpm-01.na.infn.it:1094/#g'")}'    

