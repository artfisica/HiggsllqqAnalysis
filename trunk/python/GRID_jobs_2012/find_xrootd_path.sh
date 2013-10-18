#!/bin/bash
## author: arturos@cern.ch
## Ceation: September 30th 2013
## Update:  October    7th 2013

##    INFN-t1 ds-101.cr.cnaf.infn.it oppure ds-102.cr.cnaf.infn.it (a volte ha funzionato solo ds-102.cr.cnaf.infn.it)
##    INFN-FRASCATI atlasse.lnf.infn.it
##    INFN-NAPOLI-ATLAS t2-dpm-01.na.infn.it
##    INFN-ROMA1 grid-cert-03.roma1.infn.it
##    Federated redirector per la cloud IT atlas-xrd-it.cern.ch
##    Federated redirector atlas-xrd-eu.cern.ch
##
##   REFERENCE: http://wiki.infn.it/strutture/lnf/dr/tier2/analisi/pod

###################################################################################################################################################################
FILE=$1                                                                                                                                                          ##
SITE="INFN-ROMA1_LOCALGROUPDISK" ##"INFN-NAPOLI-ATLAS_LOCALGROUPDISK"                                                                                            ##
SERVER="grid-cert-03.roma1.infn.it" ##"t2-dpm-01.na.infn.it"                                                                                                     ##
date                                                                                                                                                             ##
echo "  "                                                                                                                                                        ##
echo "################################"                                                                                                                          ##
echo "  "                                                                                                                                                        ##
                                                                                                                                                                 ##
k=1                                                                                                                                                              ##
while read files;do                                                                                                                                              ##
    while read line;do                                                                                                                                           ##
        sample=${line////}                                                                                                                                       ##
        sample=${sample//.merge.NTUP_HSG2/}                                                                                                                      ##
        echo "dq2-ls -L $SITE -fp user.$2.$sample.$3/|grep srm|awk '{print \"lcg-la \"\$1}'|awk '{system(\$0 \"|sed 's#lfn:/grid#root://$SERVER:1094/#g'\")}'"   ##
        ((k++))                                                                                                                                                  ##
    done < $files                                                                                                                                                ##   
    echo "                                     "                                                                                                                 ##
    echo " ### End of Samples = $files ###     "                                                                                                                 ##
    echo " ##################################  "                                                                                                                 ##
    echo "                                     "                                                                                                                 ##
done  < $FILE                                                                                                                                                    ##
                                                                                                                                                                 ##
echo "  "                                                                                                                                                        ##
echo "        ## Total number of Datasets use in analysis: $k ##"                                                                                                ##
echo "              ######################################     "                                                                                                 ##
echo "  "                                                                                                                                                        ##
###################################################################################################################################################################
