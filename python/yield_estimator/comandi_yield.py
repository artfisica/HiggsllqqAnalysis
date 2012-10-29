# sends yield estimations
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--step", dest="step", help="specify the step you want to start")
parser.add_option("-y", "--years", dest="years", help="specify the years you want to run on (mc12a, mc11c)")

(job_options, args) = parser.parse_args()

years = ['mc12a', 'mc11c']
flavors = ['no_constraint', 'with_constraint', 'no_constraint_no100gevcut', 'with_constraint_no100gevcut']
options_yield = ['-l', '-c -l', '', '-c']
options_replacer = ['-l', '-l', '', '']
#flavors = ['with_constraint_no100gevcut']
#options_yield = ['-c']
#options_replacer = ['']

location_data = {}
location_ZZ = {}
location_data['mc12a'] = '/afs/cern.ch/work/v/vippolit/yield_2012/data_candidates.root'
location_data['mc11c'] = '/afs/cern.ch/work/v/vippolit/yield_2011/data_candidates.root'
location_ZZ['mc12a'] = '/afs/cern.ch/work/v/vippolit/yield_2012'
location_ZZ['mc11c'] = '/afs/cern.ch/work/v/vippolit/yield_2011'

launch_dir = os.environ['PWD']
yield_estimator = 'python %s/Higgs4lepAnalysis/python/yield_estimator/yield_estimator.py' % launch_dir
signal_125_fix = 'root -b -q %s/Higgs4lepAnalysis/python/yield_estimator/Interpolate125.C' % launch_dir
signal_interpolator = '%s/HiggsZZ4lUtils/Shapes/signalInterpolation/interpolate' % launch_dir
redux_bkg_smoother = 'root -b -q %s/HiggsZZ4lUtils/Shapes/SmoothZjets/Create_smoothZxx.C+' % launch_dir
redux_bkg_smoothing_finalizer = 'root -b -q %s/HiggsZZ4lUtils/Shapes/SmoothZjets/plot.C' % launch_dir
histo_replacer = 'python %s/Higgs4lepAnalysis/python/yield_estimator/replace_histos_and_yields.py' % launch_dir

basepath = '/tmp/vippolit'

print 'INFO : called with options %s' % job_options.step

# >>> STEP 1 <<< produce MC_only
if (job_options.step == '1'):
  for year in years:
    for i, flavor in enumerate(flavors):
      directory = '%s/%s/%s' % (basepath, year, flavor)
      print 'INFO : creating MC_only folders...'
      os.system('mkdir -p %s/MC_only/workspace_histos' % (directory))
      print 'INFO : copying data...'
      os.system('cp %s %s/MC_only' % (location_data[year], directory))
      print 'INFO : launching yield estimator...'
      os.system('%s %s -d %s/MC_only -t %s &> %s/MC_only/YIELD &' % (yield_estimator, options_yield[i], directory, year, directory))

# >>> STEP 2 <<< interpolate signal shapes
elif (job_options.step == '2'):
  for year in years:
    for i, flavor in enumerate(flavors):
      directory = '%s/%s/%s' % (basepath, year, flavor)

      if (year == 'mc11c'):
        print 'WARNING : year is 2011, so we replace the 125 GeV 4e histogram with the interpolation of 120 and 130 GeV'
        os.system('cd %s/MC_only; cp %s/HiggsZZ4lUtils/Shapes/signalInterpolation/th1fmorph.C .; %s' % (directory, launch_dir, signal_125_fix));

      print 'INFO : launching signal interpolation...'

      os.system('cd %s/MC_only; ln -s __yields_tmp_outfile.root interpolationInputsZmass.root; %s &> %s/MC_only/LOG_SIGNAL_INTERPOLATION &' % (directory, signal_interpolator, directory))

# >>> STEP 3 <<< smooth reducible backgrounds
elif (job_options.step == '3'):
  for year in years:
    for i, flavor in enumerate(flavors):
      directory = '%s/%s/%s' % (basepath, year, flavor)
      os.system('mkdir -p %s/MC_only/tmp_smoothZjets' % (directory))
      print 'INFO : smoothing syst1 REDUX...'
      os.system('cd %s/MC_only/tmp_smoothZjets; cp ../reduxsyst_bkg_tree.root reduxbkg_tree.root; %s &> %s/MC_only/LOG_REDUXBKG_SMOOTHING_SYST; mv Zxx_Nominal.root Zxx_syst.root' % (directory, redux_bkg_smoother, directory))
      print 'INFO : smoothing syst2 REDUX...'
      os.system('cd %s/MC_only/tmp_smoothZjets; cp ../reduxsyst2_bkg_tree.root reduxbkg_tree.root; %s &> %s/MC_only/LOG_REDUXBKG_SMOOTHING_SYST2; mv Zxx_Nominal.root Zxx_syst2.root' % (directory, redux_bkg_smoother, directory))
      print 'INFO : smoothing nominal REDUX...'
      os.system('cd %s/MC_only/tmp_smoothZjets; cp ../reduxbkg_tree.root reduxbkg_tree.root; %s &> %s/MC_only/LOG_REDUXBKG_SMOOTHING_NOMINAL; mv Zxx_Nominal.root Zxx_nominal.root' % (directory, redux_bkg_smoother, directory))
      print 'INFO : combining to produce histograms...'
      os.system('cd %s/MC_only/tmp_smoothZjets; %s' % (directory, redux_bkg_smoothing_finalizer))

# >>> STEP 4 <<< produce the final output
elif (job_options.step == '4'):
  for year in years:
    for i, flavor in enumerate(flavors):
      directory = '%s/%s/%s' % (basepath, year, flavor)
      print 'INFO : copying smoothed ZZ...'
      os.system('cp %s/histos_*ZZSmooth_SR_2012.root %s/MC_only' % (location_ZZ[year], directory))
      print 'INFO : creating data_driven folders...'
      os.system('mkdir -p %s/data_driven/workspace_histos' % (directory))
      print 'INFO : launching histo replacer...'
      os.system('%s %s -i %s/MC_only -o %s/data_driven -t %s &> %s/data_driven/NEW_YIELDS &' % (histo_replacer, options_replacer[i], directory, directory, year, directory))
