import sys
from plots import *
from ROOT import *

#def plotinorder(plist,plot,lepchan,tagchan,norm):
#  # set options
#  lepindex = -1
#  tagindex = -1
#  if lepchan=='mumu': lepindex = 0
#  elif lepchan=='ee': lepindex = 2
#  elif lepchan=='mue': lepindex = 1
#  if tagchan=='tag': tagindex = 1
#  elif tagchan=='untag': tagindex = 0
#  elif tagchan=='incl': tagindex = 2
#  if not plot in plotlist:
#    print 'ERROR: this plot is not in plotlist'
#    sys.exit()
#
#  tmp_list = []
#  for extsample in extsamplelist:
#    tmp_list.append(plist[extsample][plot][lepindex][tagindex])
#
#  max = 0
#  maxh = 0
#  i = 0
#  opt = ''
#  leg = TLegend(.7,.7,.95,.95)
#  for h in tmp_list:
#    if norm:
#      if not h.Integral()==0: h.Scale(1./h.Integral())
#    i += 1
#    if(h.GetBinContent(h.GetMaximumBin()) >= max):
#      max = h.GetBinContent(h.GetMaximumBin())
#      maxh = i
#  canv = TCanvas()
#  canv.SetName('canvas')
#  if 'Data' in tmp_list[maxh-1].GetName(): opt = 'Pe'
#  tmp_list[maxh-1].Draw(opt)
#  i=0
#  for h in tmp_list:
#    leg.AddEntry(h,h.GetName())
#    opt = ''
#    i += 1
#    if(i!=maxh):
#      if 'Data' in h.GetName(): opt = 'Pe'
#      h.Draw(opt+'same')
#
#  leg.Draw('same')
#  canv.SaveAs(plot+'_'+lepchan+'_'+tagchan+'.root')
#
#  return
#
#def dumpyields(plist,plot,lepchan,tagchan):
#  # set options
#  lepindex = -1
#  tagindex = -1
#  if lepchan=='mumu': lepindex = 0
#  elif lepchan=='ee': lepindex = 2
#  elif lepchan=='mue': lepindex = 1
#  if tagchan=='tag': tagindex = 1
#  elif tagchan=='untag': tagindex = 0
#  elif tagchan=='incl': tagindex = 2
#  if not plot in plotlist:
#    print 'ERROR: this plot is not in plotlist'
#    sys.exit()
#
#  print ' ---- YIELDS: %s - lepchan %s - tagchan %s -----' % (plot,lepchan,tagchan)
#  for extsample in extsamplelist:
#    print '%s -- %f' % (plist[extsample][plot][lepindex][tagindex].GetName(),plist[extsample][plot][lepindex][tagindex].Integral())
#  return


## CLASS STRUCTURE

class sample:
  def __init__(self):
    self.name = ''
    self.xsec = ''
    self.entries = 0.
    self.isMC = False
    self.isBKG = False
    self.color = 0

class extendedsample:
  def __init__(self):
    self.name = ''
    self.samples = []
    self.color = 0
    self.isMC = False
    self.isBKG = False
