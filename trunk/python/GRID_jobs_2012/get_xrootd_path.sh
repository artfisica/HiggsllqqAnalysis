#!/bin/bash
## author: arturos@cern.ch
## Ceation: October 9th 2013
## Update:  October 9th 2013

##   REFERENCE: http://wiki.infn.it/strutture/lnf/dr/tier2/analisi/pod

########################################################################################################################################################################
##                                                                                                                                                                    ##
PRODUCTION=$1                                                                                                                                                         ## 
##                                                                                                                                                                    ##
rm                             xrootd_paths_$PRODUCTION.txt   path_xrootd_PoD_$PRODUCTION.txt                                                                         ##
./xrootd_HZZllqq_$PRODUCTION > xrootd_paths_$PRODUCTION.txt                                                                                                           ##
grep *.output.root             xrootd_paths_$PRODUCTION.txt > path_xrootd_PoD_$PRODUCTION.txt                                                                         ##
echo "       ------ The xrootd path for the production $PRODUCTION have been found and are saved them into the file path_xrootd_PoD_$PRODUCTION.txt ------        ";  ##
##                                                                                                                                                                    ##
########################################################################################################################################################################

