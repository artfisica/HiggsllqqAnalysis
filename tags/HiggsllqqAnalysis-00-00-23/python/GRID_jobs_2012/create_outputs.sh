#!/bin/bash
## author: arturos@cern.ch
## October 7th 2013

#############################################################################
FILE=$1                                                                    ##
                                                                           ##
date                                                                       ##
echo "  "                                                                  ##
echo "################################"                                    ##
echo "  "                                                                  ##
                                                                           ##
k=1                                                                        ##
while read files;do                                                        ##
    while read line;do                                                     ##
        sample=${line////}                                                 ##
        sample=${sample//.merge.NTUP_HSG2/}                                ##
        echo "user.$2.$sample.$3/"                                         ##
        ((k++))                                                            ##
    done < $files                                                          ##   
    echo "                       "                                         ##
    echo " ### End of Samples = $files ###     "                           ##
    echo " ##################################  "                           ##
    echo "                       "                                         ##
done  < $FILE                                                              ##
                                                                           ##
echo "  "                                                                  ##
echo "        ## Total number of Datasets produced: $k ##"                 ##
echo "              ################################     "                 ##
echo "  "                                                                  ##
#############################################################################
