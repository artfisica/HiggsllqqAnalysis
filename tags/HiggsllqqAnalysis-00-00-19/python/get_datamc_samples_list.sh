# stores in DATASET/ subfolder the data list for 2012
# author: valerio.ippolito@cern.ch

#!/bin/bash
mkdir -p DATASET

dq2-ls data12_8TeV*NTUP_HSG2*p1044_p1045/ | sort > DATASET/data12.txt # this is the 2lepton filtered sample
grep Muons DATASET/data12.txt > DATASET/data12_Muons.txt
grep Egamma DATASET/data12.txt > DATASET/data12_Egamma.txt
grep debugrec DATASET/data12.txt > DATASET/data12_debugrec.txt

dq2-ls mc12_8TeV*NTUP_HSG2*p1044/ | sort > DATASET/mc12.txt # this is the 2lepton filtered sample
