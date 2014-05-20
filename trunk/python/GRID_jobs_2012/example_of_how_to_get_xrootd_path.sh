#!/bin/bash
## author: arturos@cern.ch
## Ceation: September 30th 2013
## Update:  May  20th 2014

## User Options to Setup:
site="INFN-ROMA1_SCRATCHDISK"
dataset="user.arturos.data12_8TeV.periodD.physics_Muons.PhysCont.NTUP_2LHSG2.grp14_v01_p1649_p1650.22.0/"
rm _tmp_
dq2-ls -L $site -fp $dataset | grep srm | grep \.root >_tmp_ 
sed 's#dpm#/dpm#g'       _tmp_ > __tmp_; mv __tmp_ _tmp_
sed 's#srm://#root://#g' _tmp_ > __tmp_; mv __tmp_ _tmp_
