#!/bin/bash
## author: arturos@cern.ch
## September 30th 2013

#############################################################################
FILE=$1

if [ $4 = "yes" ]
then
   echo "        ## Tar File Creation ON... in progress now!! ##  "
   echo "  "
   echo "prun --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root\" --outDS user.$2.169039.ggH580CPS_ZZllqq.NTUP_HSG2_XXX --outputs output.root --useAthenaPackages --inDS mc12_8TeV.169039.PowhegPythia8_AU2CT10_ggH580CPS_ZZllqq.merge.NTUP_HSG2.e1622_s1581_s1586_r3658_r3549_p1344/ --nGBPerJob=10 --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --outTarBall tarball.tar --noSubmit --destSE=INFN-ROMA1_LOCALGROUPDISK"
fi

if [ $4 = "not" ]
then
   echo "        ## Tar File Creation OFF, if you need the tar, please, take it from the \"Bank of Tar\" or uncomment the line bellow ##  "
   echo "  "
   echo "## prun --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root\" --outDS user.$2.169039.ggH580CPS_ZZllqq.NTUP_HSG2_XXXX --outputs output.root --useAthenaPackages --inDS mc12_8TeV.169039.PowhegPythia8_AU2CT10_ggH580CPS_ZZllqq.merge.NTUP_HSG2.e1622_s1581_s1586_r3658_r3549_p1344/ --nGBPerJob=10 --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --outTarBall tarball.tar --noSubmit --destSE=INFN-ROMA1_LOCALGROUPDISK"
fi

echo "        ## ";date; echo "  "

echo "  "
echo "################################"
echo "  "

k=1
while read files;do
    while read line;do
        sample=${line////}
        sample=${sample//.merge.NTUP_HSG2/}
        if [ $5 = "" ]
        then
           echo "prun --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root\" --outDS user.$2.$sample.$3 --outputs output.root --useAthenaPackages --inDS $line --nGBPerJob=10 --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --inTarBall tarball.tar --destSE=INFN-ROMA1_LOCALGROUPDISK --useShortLivedReplicas"   
        else
           echo "prun --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root\" --outDS user.$2.$sample.$3 --outputs output.root --useAthenaPackages --inDS $line --nGBPerJob=10 --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --inTarBall tarball.tar --destSE=INFN-ROMA1_LOCALGROUPDISK --useShortLivedReplicas --excludedSite=$5"   
        fi
        ((k++))
    done < $files
done  < $FILE

echo "  "
echo "        ## Total number of Samples to be produced: $k ##"
echo "               ################################"
echo "  "
#############################################################################
