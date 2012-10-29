# launch needed yield estimators
#!/bin/bash

### STEP 1: produce MC_only histograms

if [ "$1" = "1" ] ; then

echo "creating MC_only folders..."
mkdir -p no_constraint/MC_only/workspace_histos
mkdir -p with_constraint/MC_only/workspace_histos
mkdir -p no_constraint_no100gevcut/MC_only/workspace_histos

echo "copying data candidates in MC_only folders..."
cp /afs/cern.ch/work/v/vippolit/yield_2011/data_candidates.root no_constraint_no100gevcut/MC_only
cp /afs/cern.ch/work/v/vippolit/yield_2011/data_candidates.root with_constraint/MC_only
cp /afs/cern.ch/work/v/vippolit/yield_2011/data_candidates.root no_constraint/MC_only

echo "running yield estimator for the MC estimate..."
python Higgs4lepAnalysis/python/yield_estimator/yield_estimator.py -d no_constraint_no100gevcut/MC_only -t mc11c &> no_constraint_no100gevcut/MC_only/YIELD &
python Higgs4lepAnalysis/python/yield_estimator/yield_estimator.py -d with_constraint/MC_only -c -l -t mc11c &> with_constraint/MC_only/YIELD &
python Higgs4lepAnalysis/python/yield_estimator/yield_estimator.py -d no_constraint/MC_only -l -t mc11c &> no_constraint/MC_only/YIELD &

fi


### STEP 2: interpolate among signals

if [ "$1" = "2" ]; then

echo "interpolating signal histograms..."
rm interpolationInputsZmass.root

ln -s no_constraint_no100gevcut/MC_only/__yields_tmp_outfile.root interpolationInputsZmass.root
HiggsZZ4lUtils/Shapes/signalInterpolation/interpolate
rm interpolationInputsZmass.root
mv interpol* no_constraint_no100gevcut/MC_only

ln -s no_constraint/MC_only/__yields_tmp_outfile.root interpolationInputsZmass.root
HiggsZZ4lUtils/Shapes/signalInterpolation/interpolate
rm interpolationInputsZmass.root
mv interpol* no_constraint/MC_only

ln -s with_constraint/MC_only/__yields_tmp_outfile.root interpolationInputsZmass.root
HiggsZZ4lUtils/Shapes/signalInterpolation/interpolate
rm interpolationInputsZmass.root
mv interpol* with_constraint/MC_only

fi

### STEP 3: produce reducible backgrounds' shapes

if [ "$1" = "3" ]; then

echo "producing smoothed REDUCIBLE background ttrees..."

basedir=$PWD
cd $basedir/no_constraint_no100gevcut/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst.root
cd $basedir/no_constraint_no100gevcut/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst2_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst2.root
cd $basedir/no_constraint_no100gevcut/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxbkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_nominal.root

cd $basedir/no_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst.root
cd $basedir/no_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst2_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst2.root
cd $basedir/no_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxbkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_nominal.root

cd $basedir/with_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst.root
cd $basedir/with_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxsyst2_bkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_syst2.root
cd $basedir/with_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; cp ../reduxbkg_tree.root reduxbkg_tree.root; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C; mv Zxx_Nominal.root Zxx_nominal.root

echo "converting to histograms..."
cd $basedir/no_constraint_no100gevcut/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/plot.C
cd $basedir/no_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/plot.C
cd $basedir/with_constraint/MC_only; mkdir -p tmp_smoothZjets; cd tmp_smoothZjets; root -q -b $TestArea/tier2/HiggsZZ4lUtils/Shapes/SmoothZjets/plot.C

cd $basedir

echo " -----> DONE !!!"

fi
### STEP 4: produce data_driven histos

if [ "$1" = "4" ]; then

echo "copying smoothed ZZ histograms..."
cp /afs/cern.ch/work/v/vippolit/yield_2011/histos_*ZZSmooth_SR_2012.root no_constraint_no100gevcut/MC_only
cp /afs/cern.ch/work/v/vippolit/yield_2011/histos_*ZZSmooth_SR_2012.root no_constraint/MC_only
cp /afs/cern.ch/work/v/vippolit/yield_2011/histos_*ZZSmooth_SR_2012.root with_constraint/MC_only

echo "creating data_driven folders..."
mkdir -p no_constraint/data_driven/workspace_histos
mkdir -p with_constraint/data_driven/workspace_histos
mkdir -p no_constraint_no100gevcut/data_driven/workspace_histos

echo "producing final results..."
python Higgs4lepAnalysis/python/yield_estimator/replace_histos_and_yields.py -i no_constraint_no100gevcut/MC_only -o no_constraint_no100gevcut/data_driven -t mc11c &> no_constraint_no100gevcut/data_driven/NEW_YIELDS &
python Higgs4lepAnalysis/python/yield_estimator/replace_histos_and_yields.py -l -i no_constraint/MC_only -o no_constraint/data_driven &> no_constraint/data_driven/NEW_YIELDS -t mc11c &
python Higgs4lepAnalysis/python/yield_estimator/replace_histos_and_yields.py -l -i with_constraint/MC_only -o with_constraint/data_driven &> with_constraint/data_driven/NEW_YIELDS -t mc11c &

echo " -----> DONE !!!"

fi
