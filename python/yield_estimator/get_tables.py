# get final yield tables from ahistos
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

# look for "kostas" in this file for those changes needed to run only in a restricted signal region

from ROOT import *

do_2mu2e_2e2mu = False
minimum_mass = 0
mass_boundary = 160
#minimum_mass = 125 # kostas
#mass_boundary = 135 # kostas

#path = '/afs/cern.ch/work/v/vippolit/kostas/latest_2011/with_constraint_no100gevcut/data_driven'
path = '/afs/cern.ch/work/v/vippolit/kostas/latest_2012/with_constraint_no100gevcut/data_driven'

masses = ['1250', '1300', '1500', '1900', '2000', '4000', '6000', '1200']

channels = ['4mu', '2mu2e', '4e', '2e2mu']

class sample:
  def __init__(self, name, histos, isSignal=False, mass=''):
    self.name = name
    self.histos = histos
    self.isSignal = isSignal
    self.mass = mass

    self.yield_tot = {}
    self.yield_low = {}
    self.yield_high = {}
    self.partial_yield_tot = {}
    self.partial_yield_low = {}
    self.partial_yield_high = {}
    for chan in channels:
      self.yield_tot[chan] = 0.
      self.yield_low[chan] = 0.
      self.yield_high[chan] = 0.
      self.partial_yield_tot[chan] = []
      self.partial_yield_low[chan] = []
      self.partial_yield_high[chan] = []
      for i in range(len(histos)):
        self.partial_yield_tot[chan].append(0)
        self.partial_yield_low[chan].append(0)
        self.partial_yield_high[chan].append(0)
      

samples = []
samples.append(sample(r'$\tabscript{ZZ}{(*)}{}$', ['histoZZ', 'histogg2ZZ']))
samples.append(sample(r'$\tabscript{qq\rightarrowZZ}{(*)}{}$', ['histoZZ'])) # for syst :)
samples.append(sample(r'$\tabscript{gg\rightarrowZZ}{(*)}{}$', ['histogg2ZZ'])) # for syst :)
samples.append(sample(r'$Z$', ['histoZ']))
samples.append(sample(r'$Zbb$', ['histoZbb']))
samples.append(sample(r'$tt$', ['histott']))
samples.append(sample(r'$Z$, $Zb\bar{b}$, and $t\bar{t}$', ['histoZ', 'histoZbb', 'histott']))
samples.append(sample(r'Total Background', ['histoZZ', 'histogg2ZZ', 'histoZ', 'histoZbb', 'histott']))
samples.append(sample(r'Data', ['histoDATA']))
samples.append(sample(r'$m_{H}=125\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '1250'))
samples.append(sample(r'$m_{H}=130\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '1300'))
samples.append(sample(r'$m_{H}=150\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '1500'))
samples.append(sample(r'$m_{H}=190\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '1900'))
samples.append(sample(r'$m_{H}=400\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '4000'))
samples.append(sample(r'$m_{H}=600\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '6000'))
samples.append(sample(r'$m_{H}=120\:\gev$', ['higgs','higgsVBF','higgsWH', 'higgsZH'], True, '1200'))
    

output_lines = []


for i, mass in enumerate(masses):
  f = TFile.Open('%s/workspace_histos/ahistos_%s.root' % (path, mass))
  print 'VERBOSE: opened %s' % f.GetName()

  # first mass: save yields
  for sample in samples:
    if ((not sample.isSignal and i == 0) or (sample.isSignal and sample.mass in mass)):
      if (sample.isSignal and sample.mass in mass): print '%s uguale a %s' % (sample.mass, mass)
      for chanid, chan in enumerate(channels):
        for histobase_idx, histobase in enumerate(sample.histos):
          histo = '%s_%d' % (histobase, chanid)

          h = f.Get(histo)
          tot = h.Integral(0, h.GetNbinsX()+1)
          low = h.Integral(h.FindBin(minimum_mass), h.FindBin(mass_boundary) - 1) # kostas
          high = h.Integral(h.FindBin(mass_boundary), h.GetNbinsX()+1)
        
          sample.yield_tot[chan] += tot
          sample.yield_low[chan] += low
          sample.yield_high[chan] += high
          sample.partial_yield_tot[chan][histobase_idx] += tot
          sample.partial_yield_low[chan][histobase_idx] += low
          sample.partial_yield_high[chan][histobase_idx] += high

          print 'INFO: retrieving %s (%lf / %lf / %lf = %lf)' % (histo, low, high, tot, low+high)

  f.Close()


for sample in samples:
  if (sample.isSignal):
    if (do_2mu2e_2e2mu):
      output_lines.append(r'%s & \multicolumn{2}{c}{%.2lf} & \multicolumn{2}{c}{%.2lf} & \multicolumn{2}{c}{%.2lf}& \multicolumn{2}{c}{%.2lf}\\' % (sample.name, sample.yield_low['4mu']+sample.yield_high['4mu'], sample.yield_low['2mu2e']+sample.yield_high['2mu2e'], sample.yield_low['2e2mu']+sample.yield_high['2e2mu'], sample.yield_low['4e']+sample.yield_high['4e']))
    else:
      output_lines.append(r'%s & \multicolumn{2}{c}{%.2lf} & \multicolumn{2}{c}{%.2lf}& \multicolumn{2}{c}{%.2lf}\\' % (sample.name, sample.yield_low['4mu']+sample.yield_high['4mu'], sample.yield_low['2mu2e']+sample.yield_high['2mu2e']+sample.yield_low['2e2mu']+sample.yield_high['2e2mu'], sample.yield_low['4e']+sample.yield_high['4e']))
  else:
    if (do_2mu2e_2e2mu):
      output_lines.append(r'%s & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf\\' % (sample.name, sample.yield_low['4mu'], sample.yield_high['4mu'], sample.yield_low['2mu2e'], sample.yield_high['2mu2e'], sample.yield_low['2e2mu'], sample.yield_high['2e2mu'], sample.yield_low['4e'], sample.yield_high['4e']))
    else:
      output_lines.append(r'%s & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf\\' % (sample.name, sample.yield_low['4mu'], sample.yield_high['4mu'], sample.yield_low['2mu2e'] + sample.yield_low['2e2mu'], sample.yield_high['2mu2e'] + sample.yield_high['2e2mu'], sample.yield_low['4e'], sample.yield_high['4e']))

print ''
print ''
print ''

for line in output_lines:
  print line

#     $\tabscript{ZZ}{(*)}{}$             & 2.47  & 19.58  & 1.78 & 29.57 & 1.25 & 12.60 \\
#     $Z$                                 & 0.00  & 0.00  & 0.59 & 0.74 & 0.51 & 4.27 \\ 
#     $Zbb$                               & 0.13  & 0.08  & 0.24 & 0.19 & 0.18 & 0.06 \\
#     $tt$                                & -0.03  & 0.03  & 0.06 & 0.00 & 0.03 & 0.06 \\\\hline
#     $Z$, $Zb\bar{b}$, and $t\bar{t}$    & 0.10  & 0.11  & 0.88 & 0.93 & 0.72 & 4.38 \\\\hline
#     Total Background                    & 2.57  & 19.69  & 2.67 & 30.51 & 1.96 & 16.98 \\\\hline 
#     Data & 2 & 24 & 2 & 52 & 3 & 19\\
#     \\hline 
#      $m_{H}=125\\:\\gev$ &\\multicolumn{2}{c}{0.98} &\\multicolumn{2}{c}{1.20}&\\multicolumn{2}{c}{0.55}\\
#      $m_{H}=130\\:\\gev$ &\\multicolumn{2}{c}{1.52} &\\multicolumn{2}{c}{1.94}&\\multicolumn{2}{c}{0.88}\\
#      $m_{H}=150\\:\\gev$ &\\multicolumn{2}{c}{3.20} &\\multicolumn{2}{c}{4.18}&\\multicolumn{2}{c}{1.92}\\
#      $m_{H}=190\\:\\gev$ &\\multicolumn{2}{c}{5.85} &\\multicolumn{2}{c}{8.89}&\\multicolumn{2}{c}{3.79}\\
#      $m_{H}=200\\:\\gev$ &\\multicolumn{2}{c}{6.36} &\\multicolumn{2}{c}{10.46}&\\multicolumn{2}{c}{4.63} \\
#      $m_{H}=400\\:\\gev$ &\\multicolumn{2}{c}{2.80} &\\multicolumn{2}{c}{4.77}&\\multicolumn{2}{c}{2.08} \\
#      $m_{H}=600\\:\\gev$ &\\multicolumn{2}{c}{0.56} &\\multicolumn{2}{c}{1.01}&\\multicolumn{2}{c}{0.44} \\
#     \\hline\\hline
#     \\hline\\hline
#     \\end{tabular}
#     \\end{center}
#     \\end{table*}
#
#%%% done %%%
#''
