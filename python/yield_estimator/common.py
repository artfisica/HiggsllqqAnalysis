from ROOT import *

def rebin_and_unbin(the_histo, n_times):
# takes an histograms and interpolates it via rebinning
  new_histo_name = (the_histo.GetName() + '_rebinned')

  result = TH1F(new_histo_name, new_histo_name, the_histo.GetNbinsX(), the_histo.GetXaxis().GetXmin(), the_histo.GetXaxis().GetXmax())

  tmp_clone = the_histo.Clone(the_histo.GetName() + '_tmp')

  nbins = tmp_clone.GetNbinsX()

  tmp_clone.Rebin(n_times)

  for i in range(0, nbins+2):
    the_point = result.GetBinCenter(i)
    result.SetBinContent(i, tmp_clone.Interpolate(the_point))

  tmp_clone.Delete()

  return result

def remove_bins_below(histo, minX):
# empties the bins below a certain values
  bins_to_empty = histo.FindBin(minX)

  for i in range(0, bins_to_empty):
    histo.SetBinContent(i, 0)

class sample:
  def __init__(self, name):
    self.generated = 0.
    self.xsec = 0.
    self.name = name
    self.histo = {}
    self.final = {}

  def __str__(self):
    lines = []
    lines.append('gen = %lf xsec = %lf\n' % (self.generated, self.xsec))
    for key in self.histo.keys():
       lines.append('[%s] int = %lf rms = %lf\n' % (key, self.histo[key].Integral(0, self.histo[key].GetNbinsX()+1), self.histo[key].GetRMS())) 

    return ''.join(lines)
