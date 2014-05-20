dq2-ls -L INFN-ROMA1_LOCALGROUPDISK -fp  user.arturos.mc12_8TeV.117293.AlpgenJimmy_AUET2CTEQ6L1_WcNp0.merge.NTUP_HSG2.e1601_s1499_s1504_r3658_r3549_p1344.13.1/|grep srm |grep \.root >_tmp_ 
sed 's#dpm#/dpm#g' _tmp_ >__tmp_; mv __tmp_ _tmp_
sed 's#srm://#root://#g' _tmp_ >__tmp_; mv __tmp_ _tmp_
