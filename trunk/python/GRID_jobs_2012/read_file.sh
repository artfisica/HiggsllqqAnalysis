#!/bin/bash
## author: arturos@cern.ch
## Created: September 30th 2013
## Updated: November  20th 2014
#############################################################################
FILE=$1
SYST=$6
PRODUCTION=$3.$6
FilesperJob="1" ## Reempacing the option: --nGBPerJob=$GbperJob
SitesToSave="INFN-ROMA1_LOCALGROUPDISK,INFN-NAPOLI-ATLAS_LOCALGROUPDISK"
CMT="--cmtConfig x86_64-slc6-gcc47-opt"
ROOTVer="--rootVer 5.34.14"
skipScout="--skipScout"
Memory="--memory=25600"

if [ $4 = "yes" ]
then
   echo "        ## Tar File Creation ON... in progress now!! ##  "
   echo "  "
   echo "prun $CMT --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root --MV1c --GSC\" --outDS user.$2.169039.ggH580CPS_ZZllqq.NTUP_HSG2_XXX --outputs output.root --inDS mc12_8TeV.169039.PowhegPythia8_AU2CT10_ggH580CPS_ZZllqq.merge.NTUP_HSG2.e1622_s1581_s1586_r3658_r3549_p1649/ --nFilesPerJob=$FilesperJob --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --outTarBall tarball.tar --noSubmit --noCompile $ROOTVer --destSE=$SitesToSave $skipScout $Memory"
fi

if [ $4 = "not" ]
then
   echo "        ## Tar File Creation OFF, if you need the tar, please, take it from the \"Bank of Tar\" or uncomment the line bellow ##  "
   echo "  "
   echo "## prun $CMT --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root --MV1c --GSC\" --outDS user.$2.169039.ggH580CPS_ZZllqq.NTUP_HSG2_XXX --outputs output.root --inDS mc12_8TeV.169039.PowhegPythia8_AU2CT10_ggH580CPS_ZZllqq.merge.NTUP_HSG2.e1622_s1581_s1586_r3658_r3549_p1649/ --nFilesPerJob=$FilesperJob --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --outTarBall tarball.tar --noSubmit --noCompile $ROOTVer --destSE=$SitesToSave $skipScout $Memory"
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
        if [ $5 = "" ] ##        if [ $5 = "" ]
        then
           echo "prun $CMT --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root --MV1c --GSC\" --outDS user.$2.$sample.$PRODUCTION --outputs output.root --inDS $line --nFilesPerJob=$FilesperJob --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --inTarBall tarball.tar --noCompile $ROOTVer --destSE=$SitesToSave --useShortLivedReplicas $skipScout $Memory"
        else
           echo "prun $CMT --exec \"echo %IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; ./HiggsllqqAnalysis/bin/test --analysis rel_17_2 --useTopoIso --input input.txt --output output.root --MV1c --GSC\" --outDS user.$2.$sample.$PRODUCTION --outputs output.root --inDS $line --nFilesPerJob=$FilesperJob --useRootCore --extFile=*/*/*.root*,*/*/*/*.root*,*/*/*/*/*.root*,*/*/*/*/*/*.root*, --inTarBall tarball.tar --noCompile $ROOTVer --destSE=$SitesToSave --useShortLivedReplicas --excludedSite=$5 $skipScout $Memory"
        fi
        ((k++))
    done < $files
done  < $FILE

echo "  "
echo "        ## Total number of Samples to be produced: $k ##"
echo "               ################################"
echo "  "
#############################################################################
