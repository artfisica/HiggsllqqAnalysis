# takes the output of yield_estimator.py and
# - replaces ZZ histos, after rebinning to [0, 1000], with histos after smoothing
# - replaces Z, Zbb and tt yields with histos after smoothing, applying data-driven estimates
# - adds interpolations for Higgs mass points (to be done)
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-t", "--tag", dest="tag", help="specify the year tag (mc11c for 2011, mc12a for 2012)", default="mc12a")
parser.add_option("-i", "--input-directory", dest="input_directory", help="expect input from directory DIR", metavar="DIR", default=".")
parser.add_option("-o", "--output-directory", dest="output_directory", help="expect output from directory DIR", metavar="DIR", default="WORKSPACES/")
parser.add_option("-l", "--lower-cut", action="store_true", default=False, help="apply lower mass cut (100 GeV)", dest="do_apply_100GeV_cut")

(options, args) = parser.parse_args()

tag = options.tag
input_dir = options.input_directory
output_dir = options.output_directory
do_apply_100GeV_cut = options.do_apply_100GeV_cut

print ''
print ''
print 'CONFIGURATION:'
print '   - tag = %s' % str(tag)
print '   - input_dir = %s' % str(input_dir)
print '   - output_dir = %s' % str(output_dir)
print '   - do_apply_100GeV_cut = %s' % str(do_apply_100GeV_cut)
print ''
print ''
print ''
print ''



def remove_bins_below(histo, minX):
# empties all bins of an histogram up to x=minX
  bins_to_empty = histo.FindBin(minX)

  for i in range(0, bins_to_empty):
    histo.SetBinContent(i, 0)


channels = ['4mu', '2mu2e', '2e2mu', '4e']

# new yields
datadriven_yield = {}
for chan in channels:
  datadriven_yield[chan] = {}

if (tag == 'mc11c'):
  datadriven_yield['4mu']['histoZ'] = 0.
  datadriven_yield['4mu']['histoZbb'] = 0.250
  datadriven_yield['4mu']['histott'] = 0.022
  datadriven_yield['2e2mu']['histoZ'] = 0.
  datadriven_yield['2e2mu']['histoZbb'] = 0.201
  datadriven_yield['2e2mu']['histott'] = 0.020
  datadriven_yield['2mu2e']['histoZ'] = 2.6
  datadriven_yield['2mu2e']['histoZbb'] = 0.
  datadriven_yield['2mu2e']['histott'] = 0.
  datadriven_yield['4e']['histoZ'] = 3.1
  datadriven_yield['4e']['histoZbb'] = 0.
  datadriven_yield['4e']['histott'] = 0.
elif (tag == 'mc12a'):
  datadriven_yield['4mu']['histoZ'] = 0.
  datadriven_yield['4mu']['histoZbb'] = 0.506
  datadriven_yield['4mu']['histott'] = 0.044
  datadriven_yield['2e2mu']['histoZ'] = 0.
  datadriven_yield['2e2mu']['histoZbb'] = 0.405
  datadriven_yield['2e2mu']['histott'] = 0.040
  datadriven_yield['2mu2e']['histoZ'] = 4.9
  datadriven_yield['2mu2e']['histoZbb'] = 0.
  datadriven_yield['2mu2e']['histott'] = 0.
  datadriven_yield['4e']['histoZ'] = 3.9
  datadriven_yield['4e']['histoZbb'] = 0.
  datadriven_yield['4e']['histott'] = 0.

# copy normalizations for shape systematics and m12, m34 histos as well
for chan in channels:
  for histoname in datadriven_yield[chan].keys():
    datadriven_yield[chan]['%s_lo' % histoname] = datadriven_yield[chan][histoname]
    datadriven_yield[chan]['%s_hi' % histoname] = datadriven_yield[chan][histoname]
    datadriven_yield[chan]['%sm12' % histoname] = datadriven_yield[chan][histoname]
    datadriven_yield[chan]['%sm34' % histoname] = datadriven_yield[chan][histoname]

chanid_name = ['4mu', '2mu2e', '4e', '2e2mu'] # maps _i to chan name

from ROOT import *
import glob

filename_smoothed_qqZZ = 'histos_qqZZSmooth_SR_2012.root'  # smoothed qq->ZZ histograms + systematics
filename_smoothed_ggZZ = 'histos_ggZZSmooth_SR_2012.root'  # smoothed gg->ZZ histograms + systematics
filename_signal_finest = '__yields_tmp_outfile.root'       # non-interpolated signal histograms (useless)
filename_interpol_sig  = 'interpol%d.root'                 # interpolated signal histograms
filename_smoothed_redux = 'tmp_smoothZjets/reducible.root' # interpolated signal histograms
filelist_ahistos_MC    = glob.glob(input_dir + '/workspace_histos/ahisto*root') # MC_only ahistos (this is the only list of files with absolute path)

list_histos_irred = ['histoZZ', 'histogg2ZZ',                     # on each row, first must be qq, then gg (so that even->qq, odd->gg)
                     'hZZ_f_shape_qcd_p', 'hgg2ZZ_f_shape_qcd_p',
                     'hZZ_f_shape_qcd_m', 'hgg2ZZ_f_shape_qcd_m',
                     'hZZ_f_shape_as_p', 'hgg2ZZ_f_shape_as_p',
                     'hZZ_f_shape_as_m', 'hgg2ZZ_f_shape_as_m',
                     'hqqZZ_E_EFF_up', 'hggZZ_E_EFF_up',
                     'hqqZZ_E_EFF_down', 'hggZZ_E_EFF_down',
                     'hqqZZ_M_EFF_up', 'hggZZ_M_EFF_up',
                     'hqqZZ_M_EFF_down', 'hggZZ_M_EFF_down',
                     'hqqZZ_E_TRIG_up', 'hggZZ_E_TRIG_up',
                     'hqqZZ_E_TRIG_down', 'hggZZ_E_TRIG_down',
                     'hqqZZ_M_TRIG_up', 'hggZZ_M_TRIG_up',
                     'hqqZZ_M_TRIG_down', 'hggZZ_M_TRIG_down',
                     'hqqZZ_E_CUTS_up', 'hggZZ_E_CUTS_up',
                     'hqqZZ_E_CUTS_down', 'hggZZ_E_CUTS_down',
                    #'hqqZZ_M_CUTS_up', 'hggZZ_M_CUTS_up',
                    #'hqqZZ_M_CUTS_down', 'hggZZ_M_CUTS_down',
		    ]
		     
list_histos_red_minimal = ['histott', 'histoZbb', 'histoZ'] # reducible backgrounds themselves
list_histos_red = ['histott', 'histoZbb', 'histoZ', # same as above plus systematics
                   'histott_lo', 'histoZbb_lo', 'histoZ_lo',
                   'histott_hi', 'histoZbb_hi', 'histoZ_hi']
list_histos_signal = ['higgs', 'higgsVBF', 'higgsWH', 'higgsZH', # signal histograms (plus energy scale systematics)
                      'higgs_ES_lo', 'higgsVBF_ES_lo', 'higgsWH_ES_lo', 'higgsZH_ES_lo',
                      'higgs_ES_hi', 'higgsVBF_ES_hi', 'higgsWH_ES_hi', 'higgsZH_ES_hi']
list_histos_data = ['histoDATA'] # data histograms
list_of_list_histos = [list_histos_irred, list_histos_red, list_histos_signal, list_histos_data]

# add m12 and m34
for this_list in list_of_list_histos:
  this_list_original = [x for x in this_list]

  for this_element in this_list_original:
    if (this_element[-3:] != '_lo' and this_element[-3:] != '_hi'): # add m12 and m34 only to non-_lo,_hi systematics
      this_list.append('%sm12' % this_element)
      this_list.append('%sm34' % this_element)
    

list_histos = [x for y in list_of_list_histos for x in y] # all the histograms

print 'INFO : following histograms will be retrieved: %s' % str(list_histos)

### CONTAINERS

# histograms which are replaced (for bkg smoothing and signal interpolation reasons)
new_histos_ZZ = {}
new_histos_redux = {}
new_histos_signal = {}

# cross-check: investigate how the signal interpolation changes those normalizations we already have
max_discrepancy = -9999
max_discrepancy_histo = 'none'

### FUNCTIONS

def set_normalization(histo, norm):
# normalizes an histogram, avoiding 0-integral ones
  previous = histo.Integral(0, histo.GetNbinsX()+1)

  if (previous != 0):
    histo.Scale(norm/previous)

  return


def get_clone_extended(histo, new_name):
# clones an histogram, possibly extending its low range from 100 to 0 GeV
  if (histo.GetXaxis().GetXmin() == 0.):
    return histo.Clone(new_name)
  else:
    print 'WARNING : %s has lowest X less than 100 GeV, extending it...' % new_name 

  binwidth = histo.GetBinWidth(1)
  n_bins_to_100 = int(100./binwidth)

  result = TH1F(new_name, '', n_bins_to_100 + histo.GetNbinsX(), 0, 1000)

  for i in range(result.GetNbinsX()):
    if (i < n_bins_to_100):
      result.SetBinContent(i, 0.)
    else:
      result.SetBinContent(i, histo.GetBinContent(i - n_bins_to_100))

  return result



### SEQUENCE

# open the reducible backgrounds' files, retrieve the smoothed histograms and save them (finest binning, widen x axis range to [0, 1000] if appropriate)
for filename in [filename_smoothed_redux]:
  f_smooth = TFile.Open('%s/%s' % (input_dir, filename))

  # retrieve histos
  h_nominal = f_smooth.Get('irredShapeNom')
  h_syst1 = f_smooth.Get('irredShapesyst1')
  h_syst2 = f_smooth.Get('irredShapesyst2')

  # save them
  for chanid in range(0, 4):
    for histo in list_histos_red_minimal:
      new_name_nominal = 'finest_%s_%d' % (histo, chanid)
      new_name_syst1 = 'finest_%s_lo_%d' % (histo, chanid)
      new_name_syst2 = 'finest_%s_hi_%d' % (histo, chanid)

     #print 'INFO: %s %s and %s being retrieved from smoothed reducible bkg histo file' % (new_name_nominal, new_name_syst1, new_name_syst2)

      new_histos_redux[new_name_nominal] = get_clone_extended(h_nominal, new_name_nominal)
      new_histos_redux[new_name_syst1]   = get_clone_extended(h_syst1, new_name_syst1)
      new_histos_redux[new_name_syst2]   = get_clone_extended(h_syst2, new_name_syst2)

      new_histos_redux[new_name_nominal].SetDirectory(0)
      new_histos_redux[new_name_syst1].SetDirectory(0)
      new_histos_redux[new_name_syst2].SetDirectory(0)

  f_smooth.Close()

print 'INFO : retrieved from REDUX file this list of smoothed histos: %s' % str(new_histos_redux.keys())

# open the gg2ZZ, qqZZ files, retrieve the smoothed histograms and save them (finest binning, widen x axis range to [0, 1000]
for i, filename in enumerate([filename_smoothed_qqZZ, filename_smoothed_ggZZ]):
  f_smooth = TFile.Open('%s/%s' % (input_dir, filename))

  for chanid in range(0, 4):
    # retrieve histos 
    h_smooth = f_smooth.Get('histoZZSmooth_%d' % chanid)          # histograms with a fixed name (same for qq and gg)
    h_syst_qcd_p = f_smooth.Get('hZZ_f_shape_qcd_p_%d' % chanid)
    h_syst_qcd_m = f_smooth.Get('hZZ_f_shape_qcd_m_%d' % chanid)
    h_syst_as_p = f_smooth.Get('hZZ_f_shape_as_p_%d' % chanid)
    h_syst_as_m = f_smooth.Get('hZZ_f_shape_as_m_%d' % chanid)

    prefix = ['hqqZZ', 'hggZZ'][i]
    h_syst_E_EFF_up = f_smooth.Get('%s_E_EFF_up_%d' % (prefix, chanid)) # histograms with a channel-dependent name (qq, gg)
    h_syst_E_EFF_down = f_smooth.Get('%s_E_EFF_down_%d' % (prefix, chanid))
    h_syst_E_TRIG_up = f_smooth.Get('%s_E_TRIG_up_%d' % (prefix, chanid)) 
    h_syst_E_TRIG_down = f_smooth.Get('%s_E_TRIG_down_%d' % (prefix, chanid))
    h_syst_E_CUTS_up = f_smooth.Get('%s_E_CUTS_up_%d' % (prefix, chanid)) 
    h_syst_E_CUTS_down = f_smooth.Get('%s_E_CUTS_down_%d' % (prefix, chanid))
    h_syst_M_EFF_up = f_smooth.Get('%s_M_EFF_up_%d' % (prefix, chanid))
    h_syst_M_EFF_down = f_smooth.Get('%s_M_EFF_down_%d' % (prefix, chanid))
    h_syst_M_TRIG_up = f_smooth.Get('%s_M_TRIG_up_%d' % (prefix, chanid)) 
    h_syst_M_TRIG_down = f_smooth.Get('%s_M_TRIG_down_%d' % (prefix, chanid))
   #h_syst_M_CUTS_up = f_smooth.Get('%s_M_CUTS_up_%d' % (prefix, chanid)) 
   #h_syst_M_CUTS_down = f_smooth.Get('%s_M_CUTS_down_%d' % (prefix, chanid))

    # their new name (basically prepends finest_ to the correct histo name, to remind it's with 0.5 GeV binning)
    new_name = ['finest_histoZZ_%d', 'finest_histogg2ZZ_%d'][i] % chanid
    new_name_syst_qcd_p = ['finest_hZZ_f_shape_qcd_p_%d', 'finest_hgg2ZZ_f_shape_qcd_p_%d'][i] % chanid
    new_name_syst_qcd_m = ['finest_hZZ_f_shape_qcd_m_%d', 'finest_hgg2ZZ_f_shape_qcd_m_%d'][i] % chanid
    new_name_syst_as_p = ['finest_hZZ_f_shape_as_p_%d', 'finest_hgg2ZZ_f_shape_as_p_%d'][i] % chanid
    new_name_syst_as_m = ['finest_hZZ_f_shape_as_m_%d', 'finest_hgg2ZZ_f_shape_as_m_%d'][i] % chanid
    new_name_syst_E_EFF_up = 'finest_%s' % (h_syst_E_EFF_up.GetName())
    new_name_syst_E_EFF_down = 'finest_%s' % (h_syst_E_EFF_down.GetName())
    new_name_syst_E_TRIG_up = 'finest_%s' % (h_syst_E_TRIG_up.GetName())
    new_name_syst_E_TRIG_down = 'finest_%s' % (h_syst_E_TRIG_down.GetName())
    new_name_syst_E_CUTS_up = 'finest_%s' % (h_syst_E_CUTS_up.GetName())
    new_name_syst_E_CUTS_down = 'finest_%s' % (h_syst_E_CUTS_down.GetName())
    new_name_syst_M_EFF_up = 'finest_%s' % (h_syst_M_EFF_up.GetName())
    new_name_syst_M_EFF_down = 'finest_%s' % (h_syst_M_EFF_down.GetName())
    new_name_syst_M_TRIG_up = 'finest_%s' % (h_syst_M_TRIG_up.GetName())
    new_name_syst_M_TRIG_down = 'finest_%s' % (h_syst_M_TRIG_down.GetName())
   #new_name_syst_M_CUTS_up = 'finest_%s' % (h_syst_M_CUTS_up.GetName())
   #new_name_syst_M_CUTS_down = 'finest_%s' % (h_syst_M_CUTS_down.GetName())

    # create the clones (possibly extending the low mass range)
    new_histos_ZZ[new_name] = get_clone_extended(h_smooth, new_name)
    new_histos_ZZ[new_name].SetDirectory(0)

    new_histos_ZZ[new_name_syst_qcd_p] = get_clone_extended(h_syst_qcd_p, new_name_syst_qcd_p)
    new_histos_ZZ[new_name_syst_qcd_p].SetDirectory(0)
    new_histos_ZZ[new_name_syst_qcd_m] = get_clone_extended(h_syst_qcd_m, new_name_syst_qcd_m)
    new_histos_ZZ[new_name_syst_qcd_m].SetDirectory(0)

    new_histos_ZZ[new_name_syst_as_p] = get_clone_extended(h_syst_as_p, new_name_syst_as_p)
    new_histos_ZZ[new_name_syst_as_p].SetDirectory(0)
    new_histos_ZZ[new_name_syst_as_m] = get_clone_extended(h_syst_as_m, new_name_syst_as_m)
    new_histos_ZZ[new_name_syst_as_m].SetDirectory(0)

    new_histos_ZZ[new_name_syst_E_EFF_up] = get_clone_extended(h_syst_E_EFF_up, new_name_syst_E_EFF_up)
    new_histos_ZZ[new_name_syst_E_EFF_up].SetDirectory(0)
    new_histos_ZZ[new_name_syst_E_EFF_down] = get_clone_extended(h_syst_E_EFF_down, new_name_syst_E_EFF_down)
    new_histos_ZZ[new_name_syst_E_EFF_down].SetDirectory(0)
    new_histos_ZZ[new_name_syst_E_TRIG_up] = get_clone_extended(h_syst_E_TRIG_up, new_name_syst_E_TRIG_up)
    new_histos_ZZ[new_name_syst_E_TRIG_up].SetDirectory(0)
    new_histos_ZZ[new_name_syst_E_TRIG_down] = get_clone_extended(h_syst_E_TRIG_down, new_name_syst_E_TRIG_down)
    new_histos_ZZ[new_name_syst_E_TRIG_down].SetDirectory(0)
    new_histos_ZZ[new_name_syst_E_CUTS_up] = get_clone_extended(h_syst_E_CUTS_up, new_name_syst_E_CUTS_up)
    new_histos_ZZ[new_name_syst_E_CUTS_up].SetDirectory(0)
    new_histos_ZZ[new_name_syst_E_CUTS_down] = get_clone_extended(h_syst_E_CUTS_down, new_name_syst_E_CUTS_down)
    new_histos_ZZ[new_name_syst_E_CUTS_down].SetDirectory(0)
    new_histos_ZZ[new_name_syst_M_EFF_up] = get_clone_extended(h_syst_M_EFF_up, new_name_syst_M_EFF_up)
    new_histos_ZZ[new_name_syst_M_EFF_up].SetDirectory(0)
    new_histos_ZZ[new_name_syst_M_EFF_down] = get_clone_extended(h_syst_M_EFF_down, new_name_syst_M_EFF_down)
    new_histos_ZZ[new_name_syst_M_EFF_down].SetDirectory(0)
    new_histos_ZZ[new_name_syst_M_TRIG_up] = get_clone_extended(h_syst_M_TRIG_up, new_name_syst_M_TRIG_up)
    new_histos_ZZ[new_name_syst_M_TRIG_up].SetDirectory(0)
    new_histos_ZZ[new_name_syst_M_TRIG_down] = get_clone_extended(h_syst_M_TRIG_down, new_name_syst_M_TRIG_down)
    new_histos_ZZ[new_name_syst_M_TRIG_down].SetDirectory(0)
   #new_histos_ZZ[new_name_syst_M_CUTS_up] = get_clone_extended(h_syst_M_CUTS_up, new_name_syst_M_CUTS_up)
   #new_histos_ZZ[new_name_syst_M_CUTS_up].SetDirectory(0)
   #new_histos_ZZ[new_name_syst_M_CUTS_down] = get_clone_extended(h_syst_M_CUTS_down, new_name_syst_M_CUTS_down)
   #new_histos_ZZ[new_name_syst_M_CUTS_down].SetDirectory(0)

  f_smooth.Close()

print 'INFO : retrieved from ZZ file this list of smoothed histos: %s' % str(new_histos_ZZ.keys())

# open the interpolation files for signal, retrieve the histograms and save them (finest binning)
for chanid in range(0, 4):
  f_interpol = TFile.Open('%s/%s' % (input_dir, filename_interpol_sig % chanid))

  keys_in_file = f_interpol.GetListOfKeys()

  for k in keys_in_file: 
    plot_name = k.GetName()

    h_interpol = f_interpol.Get(plot_name)

    new_histos_signal[plot_name] = h_interpol.Clone()
    new_histos_signal[plot_name].SetDirectory(0)

  f_interpol.Close()

print 'INFO : retrieved from SIGNAL file this list of interpolated histos: %s' % str(new_histos_signal.keys())

print ''
print ''
print ''


# loop over each input (and so output) ahistos file (for normalization and additional histograms like m12, m34)
for ahisto in filelist_ahistos_MC:
  f_ahisto = TFile.Open(ahisto)
  print 'VERBOSE : opening %s' % ahisto

  # this dictionary will contain all the histograms in this file
  list_histos_in_file = {} # dictionary key is the final name of the histo (i.e. histoZ_0 not finest_histoZ_0)

  # retrieve the original (i.e. MC_only) interesting histograms in this ahisto and clone them
  for histo in list_histos:
    for chanid in range(0, 4):
      histoname = '%s_%d' % (histo, chanid)
      tmp_histo = f_ahisto.Get(histoname)
      if (tmp_histo):
        # clone the histogram with a different name
        list_histos_in_file[histoname] = tmp_histo.Clone('clone_%s' % tmp_histo.GetName())
        list_histos_in_file[histoname].SetDirectory(0)
      else:
        print 'WARNING : unable to retrieve the expected histogram %s !!!' % histoname

  # close the input file
  f_ahisto.Close()

  # retrieve the binning needed for this mass hypothesis
  binning = list_histos_in_file['histoZZ_0'].GetBinWidth(1)

  # retrieve the yields for ZZ
  qqZZ_yield = {}
  ggZZ_yield = {}

  for chanid in range(0, 4):
    name_qqZZ = '%s_%d' % (list_histos_irred[0], chanid) # remember the convention of list_histos_irred filling: even position = qq, odd = gg
    name_ggZZ = '%s_%d' % (list_histos_irred[1], chanid)

    qqZZ_yield[chanid] = list_histos_in_file[name_qqZZ].Integral(0, list_histos_in_file[name_qqZZ].GetNbinsX()+1)
    ggZZ_yield[chanid] = list_histos_in_file[name_ggZZ].Integral(0, list_histos_in_file[name_ggZZ].GetNbinsX()+1)

    # delete the original (non-smoothed) ZZ histograms
    list_histos_in_file[name_qqZZ].Delete()
    list_histos_in_file[name_ggZZ].Delete()

    # but save m12 and m34 ZZ histograms with the proper name (they are not smoothed)
    name_qqZZ_m12 = '%sm12_%d' % (list_histos_irred[0], chanid)
    name_ggZZ_m12 = '%sm12_%d' % (list_histos_irred[1], chanid)
    name_qqZZ_m34 = '%sm34_%d' % (list_histos_irred[0], chanid)
    name_ggZZ_m34 = '%sm34_%d' % (list_histos_irred[1], chanid)

    list_histos_in_file[name_qqZZ_m12].SetName(name_qqZZ_m12)
    list_histos_in_file[name_ggZZ_m12].SetName(name_ggZZ_m12)
    list_histos_in_file[name_qqZZ_m34].SetName(name_qqZZ_m34)
    list_histos_in_file[name_ggZZ_m34].SetName(name_ggZZ_m34)


  # open the output file in output_dir
  f_ahisto_new = TFile('%s/workspace_histos/%s' % (output_dir, ahisto.split('/')[-1]), 'RECREATE')

  # creates the proper histogram for ZZ (clones the smoothed ones, normalizes them and rebins to match this mass hypothesis' binwidth)
  for chanid in range(0, 4):
    for h_index, histo in enumerate(list_histos_irred): # use enumerate (again, 2k -> qq, 2k+1 -> gg)
      name_histo = '%s_%d' % (histo, chanid)

      if 'finest_%s' % name_histo not in new_histos_ZZ.keys(): # do not rebin or normalize the non-smoothed ZZ histograms
        print 'WARNING : %s not in new_histos_ZZ -> will not be rebinned and normalized manually!' % name_histo
	continue
      else:
        print 'INFO : normalizing %s to the ZZ expected yield' % name_histo

      # so, clone (note that the original list_histos_in_file[name_histo] was deleted)
      if (name_histo in list_histos_in_file.keys()):
        print 'ERROR : %s is among the keys of list_histos_in_file, but we should have deleted it !!! value is %s' % (name_histo, str(list_histos_in_file[name_histo]))
      list_histos_in_file[name_histo] = new_histos_ZZ['finest_%s' % name_histo].Clone(name_histo)

      # normalize
      its_area = qqZZ_yield[chanid] if (h_index % 2 == 0) else ggZZ_yield[chanid] # still the same convention of course
      set_normalization(list_histos_in_file[name_histo], its_area)

      # rebin
      list_histos_in_file[name_histo].Rebin(int(binning / list_histos_in_file[name_histo].GetBinWidth(1)))

      # empty the low-mass bins if appropriate
      if (do_apply_100GeV_cut):
        remove_bins_below(list_histos_in_file[name_histo], 100.)

  # normalizes properly the reducible backgrounds, changing the name accordingly
  for chanid in range(0, 4):
    for histo in list_histos_red: # contains also _lo and _hi of course
      name_histo = '%s_%d' % (histo, chanid)

      if 'finest_%s' % name_histo in new_histos_redux.keys():
        # replace the histo with the new one after smoothing
	print 'INFO : replacing %s with the smoothed version' % name_histo
        list_histos_in_file[name_histo].Delete()
        list_histos_in_file[name_histo] = new_histos_redux['finest_%s' % name_histo].Clone(name_histo)
      else:
        # replace just name
	print 'INFO: keeping %s at the non-smoothed version' % name_histo
        list_histos_in_file[name_histo].SetName(name_histo)

      # normalize properly
      set_normalization(list_histos_in_file[name_histo], datadriven_yield[chanid_name[chanid]][histo])
      list_histos_in_file[name_histo].Rebin(int(binning / list_histos_in_file[name_histo].GetBinWidth(1)))

      # empty the low-mass bins if appropriate
      if (do_apply_100GeV_cut):
        remove_bins_below(list_histos_in_file[name_histo], 100.)

  # act on signals
  for chanid in range(0, 4):
    for histo in list_histos_signal:
      # name of the interpolated and the final histogram
      name_histo = '%s_%d_inter_%s' % (histo, chanid, ahisto[-9:-5]) # ahisto ends by e.g. _1300.root
      name_histo_new = '%s_%d' % (histo, chanid)

      if (name_histo not in new_histos_signal.keys()):
        print 'ERROR : interpolation not found for %s, using default histo!!!' % name_histo

        list_histos_in_file[name_histo_new].SetName(name_histo_new) # no need to rebin, it's already okay
      else:
        # keep track of the discrepancies between interpolated and non-interpolated yields
	if (list_histos_in_file[name_histo_new].Integral() != 0.):
  	  this_discrepancy = abs((list_histos_in_file[name_histo_new].Integral() - new_histos_signal[name_histo].Integral())/list_histos_in_file[name_histo_new].Integral())
	else:
	  this_discrepancy = -1
	if this_discrepancy > max_discrepancy:
	  max_discrepancy = this_discrepancy
	  max_discrepancy_histo = name_histo

        # replace the histogram found in the file with the interpolated one
        list_histos_in_file[name_histo_new].Delete() # one does not need the original histo
        list_histos_in_file[name_histo_new] = new_histos_signal[name_histo] # take the interpolated
        list_histos_in_file[name_histo_new].SetName(name_histo_new) # give it the final name
        list_histos_in_file[name_histo_new].Rebin(int(binning / list_histos_in_file[name_histo_new].GetBinWidth(1))) # they must be rebinned since interpol is at 0.5 gev

	# empty the low-mass bins if appropriate
        if (do_apply_100GeV_cut):
          remove_bins_below(list_histos_in_file[name_histo_new], 100.)

  # act on data
  for chanid in range(0, 4):
    for histo in list_histos_data:
      name_histo = '%s_%d' % (histo, chanid)
      list_histos_in_file[name_histo].SetName(name_histo) # it should be useless 

      # empty the low-mass bins if appropriate
      if (do_apply_100GeV_cut):
        remove_bins_below(list_histos_in_file[name_histo], 100.)

  # set directory and save
  print 'INFO : ******* '
  print 'INFO : %s [%s]' % (ahisto, chanid_name[chanid])

  for histo in list_histos_in_file:
    list_histos_in_file[histo].SetDirectory(f_ahisto_new)

    # verbose
    print 'INFO : %s has yield %lf and binwidth %lf' % (histo, list_histos_in_file[histo].Integral(0, list_histos_in_file[histo].GetNbinsX()+1), list_histos_in_file[histo].GetBinWidth(1))

  f_ahisto_new.Write()
  f_ahisto_new.Close()



# debug
print ''
print 'WARNING : maximum discrepancy between interpolated and original (relative) signal histograms is %lf (histo=%s)' % (max_discrepancy, max_discrepancy_histo)
