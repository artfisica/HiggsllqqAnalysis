##################################################################
# python script to make plots from the
# TestSelection output ntuples
# author Francesco Lo Sterzo <francesco.lo.sterzo@cern.ch>
##################################################################

############### INFO ###########################################
#
# USAGE python {TestSelection/python/}plotterino.py 
#
# the list of exteded dataset used is in extendedsamples.py
#
# plot strucrure is a multi-dimentional dictionary
# plot[extsample][plot][mumu/mue/ee][tag/untag/inclusive] = TH1D()
###############################################################

import array
import PyCintex
import re
import sys
PyCintex.Cintex.Enable()
from ROOT import *
gROOT.SetBatch(True)
gROOT.SetStyle('Plain')

from samples import *
from extendedsamples import *
from plots import *

############### USEFUL FUNCTIONS ###################

def plotinorder(plist,plot,lepchan,tagchan,norm):
  # set options
  lepindex = -1
  tagindex = -1
  if lepchan=='MU2': lepindex = 0
  elif lepchan=='E2': lepindex = 2
  elif lepchan=='MUE': lepindex = 1
  if tagchan=='tag': tagindex = 1
  elif tagchan=='untag': tagindex = 0
  elif tagchan=='incl': tagindex = 2
  if not plot in plotlist:
    print 'ERROR: this plot is not in plotlist'
    sys.exit()

  tmp_list = []
  sig = plist[Signal_ggH][plot][lepindex][tagindex].Clone('Signal')
  for extsample in extsamplelist:
    if not 'Signal' in extsample.name:
      tmp_list.append(plist[extsample][plot][lepindex][tagindex])
    elif extsample.name=='Signal_VBFH':
      sig.Add(plist[extsample][plot][lepindex][tagindex])
      tmp_list.append(sig)

  max = 0
  maxh = 0
  i = 0
  opt = ''
  leg = TLegend(.7,.7,.95,.95)
  for h in tmp_list:
    if norm:
      if not h.Integral()==0: h.Scale(1./h.Integral())
    i += 1
    if(h.GetBinContent(h.GetMaximumBin()) >= max):
      max = h.GetBinContent(h.GetMaximumBin())
      maxh = i
  canv = TCanvas()
  canv.SetName('canvas')
  if 'Data' in tmp_list[maxh-1].GetName(): opt = 'Pe'
  tmp_list[maxh-1].Draw(opt)
  i=0
  for h in tmp_list:
    leg.AddEntry(h,h.GetName())
    opt = ''
    i += 1
    if(i!=maxh):
      if 'Data' in h.GetName(): opt = 'Pe'
      h.Draw(opt+'same')

  leg.Draw('same')
  canv.SaveAs(plot+'_'+lepchan+'_'+tagchan+'.root')
  canv.SaveAs(plot+'_'+lepchan+'_'+tagchan+'.eps')
  canv.SetLogy(1)
  canv.SaveAs(plot+'_'+lepchan+'_'+tagchan+'_log.eps')

  return








def dumpyields(plist,plot,lepchan,tagchan):
  # set options
  lepindex = -1
  tagindex = -1
  if lepchan=='MU2': lepindex = 0
  elif lepchan=='E2': lepindex = 2
  elif lepchan=='MUE': lepindex = 1
  if tagchan=='tag': tagindex = 1
  elif tagchan=='untag': tagindex = 0
  elif tagchan=='incl': tagindex = 2
  if not plot in plotlist:
    print 'ERROR: this plot is not in plotlist'
    sys.exit()

  print ' ---- YIELDS: %s - lepchan %s - tagchan %s -----' % (plot,lepchan,tagchan)
  for extsample in extsamplelist:
    print '%s -- %f' % (plist[extsample][plot][lepindex][tagindex].GetName(),plist[extsample][plot][lepindex][tagindex].Integral(0,plist[extsample][plot][lepindex][tagindex].GetNbinsX()+1))
  return










def dofractionfit(hdata,hmc,hqcd,ndata,nmc,nqcd):
	ndata = hdata.Integral()
	nmc.append(hmc.Integral())
	nqcd.append(hqcd.Integral())
	fqcd = Double(1)
	eqcd = Double(1)
	fmc = Double(1)
	emc = Double(1)
	templates = TObjArray(2)
	templates.Add(hqcd)
	templates.Add(hmc)
	thefitter = TFractionFitter(hdata,templates)
	thefitter.Constrain(0,0.,1.)
	thefitter.Constrain(1,0.,1.)
        ## thefitter.SetRangeX(8,32)  ## 20-80GeV in mll plot  (2.5GeV per bin) 
        
	thefit = thefitter.GetFitter()
	thefit.SetParameter(0,'fqcd',.0005,.0001,0.,1.)
	thefit.SetParameter(1,'fmc',.9995,.0001,0.,1.)
	status = thefitter.Fit()
	if status==0:
		thefitter.GetResult(0,fqcd,eqcd)
		thefitter.GetResult(1,fmc,emc)
	nmc.append(fmc)
	nmc.append(emc)
	nqcd.append(fqcd)
	nqcd.append(eqcd)

	print ' -- FIT RESULT: sample: %s | status = %d ' % (hdata.GetName(),status)
	print ' qcd fraction = %f +/- %f' % (nqcd[1],nqcd[2])
	print ' mc fraction  = %f +/- %f' % (nmc[1],nmc[2])
	print ' chisquare = %f/%d' % (thefitter.GetChisquare(),thefitter.GetNDF())
	## plot the residuals
	reshisto = TH1D('reshisto','data-fit',100,-15,15)
	for b in range(hdata.GetNbinsX()): 
		fitvalue = hmc.GetBinContent(b+1)*fmc + hqcd.GetBinContent(b+1)*fqcd
		res = (hdata.GetBinContent(b+1)-fitvalue)
		reshisto.Fill(res)
	reshisto.SaveAs(hdata.GetName()+'.root')

	return







def scaleQCD(plotlist,plot,lepchan,tagchan,fmc,fqcd):
  # set options
  lepindex = -1
  tagindex = -1
  if lepchan=='MU2': lepindex = 0
  elif lepchan=='E2': lepindex = 2
  elif lepchan=='MUE': lepindex = 1
  if tagchan=='tag': tagindex = 1
  elif tagchan=='untag': tagindex = 0
  elif tagchan=='incl': tagindex = 2


## MARCO'S WAY
  rqcd = fqcd
  rmc = fmc
  for extsample in extsamplelist:
    if extsample.name=='QCD':
      plotlist[extsample][plot][lepindex][tagindex].Scale(rqcd)
#    elif extsample.isMC:
#      plotlist[extsample][plot][lepindex][tagindex].Scale(rmc)

  print '%s %s %s - QCD \'SFs\': \n Rmc = %f (fmc = %f)\n Rqcd = %f (fqcd = %f)\n' % (plot,lepchan,tagchan,rmc,fmc,rqcd,fqcd)
  return plotlist








def calculateSFz(plist,plot,lepchan,tagchan):
  # set options
  lepindex = -1
  tagindex = -1
  if lepchan=='MU2': lepindex = 0
  elif lepchan=='E2': lepindex = 2
  elif lepchan=='MUE': lepindex = 1
  if tagchan=='tag': tagindex = 1
  elif tagchan=='untag': tagindex = 0
  elif tagchan=='incl': tagindex = 2

  data = plist[Data][plot][lepindex][tagindex].Integral(0,plist[Data][plot][lepindex][tagindex].GetNbinsX()+1)
  Z = plist[Zlf_Alp][plot][lepindex][tagindex].Integral(0,plist[Zlf_Alp][plot][lepindex][tagindex].GetNbinsX()+1)
  MC = 0.
  for extsample in extsamplelist:
	  if extsample.isBKG:
		  if (not extsample.name=='Zlf_Alp') and (not extsample.name=='totbkg'): MC+= plist[extsample][plot][lepindex][tagindex].Integral(0,plist[extsample][plot][lepindex][tagindex].GetNbinsX()+1)

  print ' Data = %f; MC = %f; Z = %f;    plist = %s, plot = %s, lepchan = %s, tagchan = %s \n' % (data,MC,Z,plist,plot,lepchan,tagchan)

  return (data-MC)/Z
##################################################








if len(sys.argv)<3:
  print 'ERROR: wrong arguments:\n right usage: python TestSelection/python/plotterino.py <MU2/MUE/E2> <higgs-mass>'
  sys.exit()
chan = sys.argv[1]
mass = sys.argv[2]

add_tag = ''
lumi = 13000. # [1/pb]
data_qcd_flag = ''

qcdflag = [
  '',
  '',
  '',
]

###
# "default" QCD samples are: SS, isolated (Muons) - tightPP iso-not iso (XOR), OS (Electrons)
qcdflag[0] = ' && isqcdevent==1 && lep1_charge*lep2_charge==1' ## SS, isolated muon selection
## qcdflag[0] = ' && isqcdevent==1 && lep1_charge*lep2_charge==1 && lep1_trackiso<.1 && TMath::Abs(lep1_d0)/lep1_sigd0<3.5 && lep2_trackiso<.1 && TMath::Abs(lep2_d0)/lep2_sigd0<3.5' ## SS, isolated muon selection
#qcdflag[0] = ' && isqcdevent==1 && (lep1_trackiso > .3 || lep2_trackiso > .3) && TMath::Abs(lep1_d0)/lep1_sigd0<3.5 && TMath::Abs(lep2_d0)/lep2_sigd0<3.5 && lep1_charge*lep2_charge==-1 ' ## iso-notiso muon selection
qcdflag[1] = ' && isqcdevent==1 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1)) && lep1_charge*lep2_charge==-1' ## emu
## qcdflag[2] = ' && isqcdevent==1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1)) && lep1_charge*lep2_charge==-1 && lep1_quality==3 && lep2_quality==3' ## tightPP iso-notiso (XOR) OS electrons
## qcdflag[2] = ' && isqcdevent==1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1)) && lep1_charge*lep2_charge==1 && lep1_quality==3 && lep2_quality==3' ## tightPP iso-notiso (XOR) SS electrons
## qcdflag[2] = ' && isqcdevent==1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1)) && lep1_charge*lep2_charge==-1 && lep1_quality==2 && lep2_quality==2' ## mediumPP iso-notiso (XOR) OS electrons
#qcdflag[2] = ' && isqcdevent==1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1)) && lep1_charge*lep2_charge==-1 && lep1_quality==1 && lep2_quality==1' ## loosePP iso-notiso (XOR) OS electrons
qcdflag[2] = ' && isqcdevent==1 && (lep1_charge*lep2_charge)==1 && (TMath::Abs(lep1_d0)/lep1_sigd0)>6.5 && (TMath::Abs(lep1_d0)/lep1_sigd0)>6.5 && ((lep1_trackiso > .3 && lep2_trackiso < .1) || (lep2_trackiso > .3 && lep1_trackiso < .1))'  ##New selection 2012



dataflag = [
  '',
  '',
  ''
]

dataflag[0] = ' && isqcdevent==0 && lep1_charge*lep2_charge==-1 && lep1_trackiso < .1 && lep2_trackiso < .1 && TMath::Abs(lep1_d0)/lep1_sigd0 < 3.5 && TMath::Abs(lep2_d0)/lep2_sigd0 < 3.5' ## mumu
dataflag[1] = ' && isqcdevent==0 && lep1_charge*lep2_charge==-1' ## mue
dataflag[2] = ' && isqcdevent==0 && lep1_charge*lep2_charge==-1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && lep1_trackiso < .1 && lep2_trackiso < .1' ## ee

mcflag = [
  '',
  '',
  '',
]

mcflag[0] = '  && lep1_charge*lep2_charge==-1 && lep1_trackiso < .1 && lep2_trackiso < .1 && TMath::Abs(lep1_d0)/lep1_sigd0 < 3.5 && TMath::Abs(lep2_d0)/lep2_sigd0 < 3.5' ##mumu
mcflag[1] = '  && lep1_charge*lep2_charge==-1' ##mue
mcflag[2] = '  && lep1_charge*lep2_charge==-1 && TMath::Abs(lep1_d0)/lep1_sigd0<6.5 && TMath::Abs(lep2_d0)/lep2_sigd0<6.5 && lep1_trackiso < .1 && lep2_trackiso < .1' ##ee



tagcondition = [
  'isTagged==0',
  'isTagged==1',
  'isTagged!=-1',
]



lepindexlist = []
if chan=='E2': lepindexlist = [2]
elif chan=='MU2': lepindexlist = [0]
elif chan=='MUE': lepindexlist = [1]

tagindexlist = range(3)





extsamplelist = [
  Zlf_Alp,
  Wlf_Alp,
  Top,
  Diboson,
  # DYlf_Alp,
  # DYhf_Alp,
  # Zhf_Alp,
  Signal_ggH,
  Signal_VBFH,
  QCD,
  TotalBkg,
  Data,
]




if mass=='120':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH120_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH110_ZZllqq)
elif mass=='125':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH125_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH125_ZZllqq)
elif mass=='130':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH130_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH130_ZZllqq)
elif mass=='135':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH135_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH135_ZZllqq)
elif mass=='140':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH140_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH140_ZZllqq)
elif mass=='145':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH145_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH145_ZZllqq)
elif mass=='150':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH150_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH150_ZZllqq)
elif mass=='155':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH155_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH155_ZZllqq)
elif mass=='160':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH160_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH160_ZZllqq)
elif mass=='165':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH165_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH165_ZZllqq)
elif mass=='170':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH170_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH170_ZZllqq)
elif mass=='175':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH175_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH175_ZZllqq)
elif mass=='180':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH180_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH180_ZZllqq)
elif mass=='400':
  Signal_ggH.samples.append(PowhegPythia8_AU2CT10_ggH400_ZZllqq)
  Signal_VBFH.samples.append(PowhegPythia8_AU2CT10_VBFH400_ZZllqq)
else:
  print 'ERROR: wrong Higgs mass value. You can choose between 120 and 180 with 5 GeV steps\n'
  sys.exit()



plot = {}
for extsample in extsamplelist:
  plot[extsample] = {}
  for p in plotlist: 
    plot[extsample][p] = {}
    for lep in lepindexlist: ## mumu/emu/ee
      plot[extsample][p][lep] = {}
      for tag in tagindexlist: ## untag/tag/inclusive
        plot[extsample][p][lep][tag] = TH1D()
print 'Plot dictionary properly initialized'




# read the input files from input.txt
input_txt = open('./local_3.txt','r')
filelist = input_txt.readlines()
input_txt.close()



nqcdtag = 0.
nqcdtemplate = 0.




for extsample in extsamplelist:
  for p in plotlist:
    for lepindex in lepindexlist:
      for tagindex in tagindexlist:
        ## WARNING: need to handle this properly
         plot[extsample][p][lepindex][tagindex] = TH1D(plotlist[plotlist.index(p)]+'_'+lepchan[lepindex]+'_'+tagchan[tagindex]+'_'+extsample.name,'',nbins[plotlist.index(p)],binmin[plotlist.index(p)],binmax[plotlist.index(p)])
      ## end of loop on tag
    ## end of loop in lep chans
  ## end of loop on plots

  for sample in extsample.samples:
    sample.entries = 1.
#    file_regexp = re.compile('.*/user.*'+sample.name+'.*'+add_tag+'.*/.*root')
    file_regexp = re.compile('.*/'+sample.name+'_'+chan+'.root')
    filename = ''
    for file in filelist:
      foundfile = file_regexp.search(file)
      if foundfile:
        tmp = file
        filename = tmp.rstrip('\n')
    ## need to insert a check against more then 1 files available !!!!
    if filename=='':
      print 'ERROR: file not found for sample %s' % sample.name
      continue
    f = TFile.Open(filename)
    print 'DEBUG: %s %s ' % (extsample.name,filename)
    cutflow = f.Get('generatedEntriesHisto') ## ('cutflow_weight')
    tmpcf = f.Get('generatedEntriesHisto')   ## ('cutflow')
    print 'DEBUG: %s %s %f %f' % (extsample.name,filename,cutflow.GetBinContent(1),tmpcf.GetBinContent(1))
    sample.entries += cutflow.GetBinContent(1)
    tree = f.Get('tree')
    for p in plotlist:
      for lepindex in lepindexlist:
        for tagindex in tagindexlist:
          h_tmp = plot[extsample][p][lepindex][tagindex].Clone('tmp')
          h_tmp.Reset('MICE')
          if extsample.name=='Data': data_qcd_flag = dataflag[lepindex]
          elif extsample.name=='QCD': data_qcd_flag = qcdflag[lepindex]
          else: data_qcd_flag = mcflag[lepindex]
          tmptagindex = tagindex
          if extsample.name=='QCD': tmptagindex=2 ## to be uncommented once the normalization is OK
          condition = '('+basecondition[plotlist.index(p)]+' && channel=='+str(lepindex)+' && '+tagcondition[tmptagindex]+data_qcd_flag+' && weight>=0)' ## *weight*ggFweight'+btagsf[plotlist.index(p)]+analysis_selection[plotlist.index(p)]
          ## tmp - before trig_SF correction
          if extsample.isMC: condition = condition + '*trig_SF'
          print 'DEBUG: %s %s %s %s ' % (extsample.name,filename,p,whattodraw[plotlist.index(p)])
          tree.Draw(whattodraw[plotlist.index(p)]+'>>tmp',condition)
          if p=='mllqq_High': print condition
          if extsample.name=='QCD':
            print 'DEBUG: plot = %s plotindex = %d lepindex = %d tagindex = %d' % (p,plotlist.index(p),lepindex,tagindex)
            print condition
            print '-------------------------------------'
          ## additional normalization factor in the tagged channel
          if extsample.name=='QCD' and p=='mll_2j_High': 
            nqcdtemplate = h_tmp.Integral(0,h_tmp.GetNbinsX()+1)
            print 'arturo = %f' % (h_tmp.Integral(0,h_tmp.GetNbinsX()+1))
            if tagindex==1:
              hhh = TH1D('hhh','',100,0,200000)
              tree.Draw(whattodraw[plotlist.index(p)]+'>>hhh','('+basecondition[plotlist.index(p)]+' && channel=='+str(lepindex)+' && isTagged==1 && weight>=0'+data_qcd_flag+')') ## *weight'+btagsf[plotlist.index(p)]+analysis_selection[plotlist.index(p)])
              nqcdtag = hhh.Integral(0,h_tmp.GetNbinsX()+1)
              print 'INFO: plot %s whattodraw %s tagindex %d' % (p,whattodraw[plotlist.index(p)],tagindex)
              print 'QCD condition: %s' % '('+basecondition[plotlist.index(p)]+' && channel=='+str(lepindex)+' && isTagged==1 && weight>=0'+data_qcd_flag+')' ## *weight'+btagsf[plotlist.index(p)]+analysis_selection[plotlist.index(p)]
              print 'ADDITIONAL QCD SF: nqcdtag = %f nqcdtemplate = %f SF = %f' % (nqcdtag,nqcdtemplate,nqcdtag/nqcdtemplate)
#          if p=='mll_2j_High': print '%s - %s - %s' % (extsample.name,sample.name,condition)
          scale = 1.
          if sample.isMC:
            if sample.entries>0: scale = sample.xsec / sample.entries * lumi
            h_tmp.Scale(scale)
          plot[extsample][p][lepindex][tagindex].Add(h_tmp)
        ## end of loop on tag
      ## end of loop in lep chans
    ## end of loop on plots
    f.Close()
  ## end of loop on samples
  ## some style here
  for p in plotlist:
    for lepindex in lepindexlist:
      for tagindex in tagindexlist:
        if extsample.isMC or extsample.isBKG:
          plot[extsample][p][lepindex][tagindex].SetLineWidth(2)
          plot[extsample][p][lepindex][tagindex].SetLineColor(extsample.color)
          if not extsample.isBKG:
            plot[extsample][p][lepindex][tagindex].SetFillStyle(3002)
            plot[extsample][p][lepindex][tagindex].SetFillColor(extsample.color)
        ## if isMC or isBKG
        else:
          plot[extsample][p][lepindex][tagindex].SetMarkerStyle(23)
        ## data
      ##
    ##
  ## 
## end of loop on extended samples








##### tagged channel rebinning ########
lepindex = lepindexlist[0]
tagindex = 1 ## tag
for extsample in extsamplelist:
  plot[extsample]['mllqq_High'][lepindex][tagindex].Rebin(2)
  plot[extsample]['mllqq_sb_High'][lepindex][tagindex].Rebin(2)
  plot[extsample]['mllqq_lsb_High'][lepindex][tagindex].Rebin(2)
  plot[extsample]['mllqq_hsb_High'][lepindex][tagindex].Rebin(2)
######################################






## Fit the QCD fraction
tagindex = 2 #inclusive
lepindex = lepindexlist[0]
whichplot = 'mll_2j_High'
tff_qcd_2j = plot[QCD][whichplot][lepindex][tagindex].Clone('qcd_2j')
tff_data_2j = plot[Data][whichplot][lepindex][tagindex].Clone('data_2j')
tff_totmc_2j = plot[QCD][whichplot][lepindex][tagindex].Clone('totmc_2j')
tff_totmc_2j.Reset('MICE')

###
for extsample in extsamplelist:
	if extsample.isMC and not extsample.name=='totbkg': tff_totmc_2j.Add(plot[extsample][whichplot][lepindex][tagindex])

###
whichplot = 'mll_High'
tff_qcd = plot[QCD][whichplot][lepindex][tagindex].Clone('qcd')
tff_data = plot[Data][whichplot][lepindex][tagindex].Clone('data')
tff_totmc = plot[QCD][whichplot][lepindex][tagindex].Clone('totmc')
tff_totmc.Reset('MICE')

###
for extsample in extsamplelist:
	if extsample.isMC and not extsample.name=='totbkg': tff_totmc.Add(plot[extsample][whichplot][lepindex][tagindex])



ndata = 0.
nmc = []
nqcd = []


ndata_2j = 0.
nmc_2j = []
nqcd_2j = []



dofractionfit(tff_data,tff_totmc,tff_qcd,ndata,nmc,nqcd)
dofractionfit(tff_data_2j,tff_totmc_2j,tff_qcd_2j,ndata_2j,nmc_2j,nqcd_2j)



n_index = 0 # integral before fit
f_index = 1 # fraction (from fit)
e_index = 2 # error on the fraction (form fit)



## plot = scaleQCD(plot,'chisquare',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mll',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
## plot = scaleQCD(plot,'mll_2j',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq',chan,'tag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral()*nqcdtag/nqcdtemplate)
## plot = scaleQCD(plot,'mllqq',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'leadleppt',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
## plot = scaleQCD(plot,'subleppt',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
## plot = scaleQCD(plot,'mllqq_sb',chan,'tag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral()*nqcdtag/nqcdtemplate)
## plot = scaleQCD(plot,'mllqq_sb',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mjj',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_lsb',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_hsb',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())

################################
## plot = scaleQCD(plot,'chisquare_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mll_Low',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
## plot = scaleQCD(plot,'ptll_Low',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
## plot = scaleQCD(plot,'mll_2j_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'leadleppt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'subleppt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'leadlepeta_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'sublepeta_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1pt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2pt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'corrJ1pt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'corrJ2pt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realZpt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'corrZpt_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realZm_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'corrZm_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'nJets_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'nbJets_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_MV1_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_jvf_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_ntrk_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## ## plot = scaleQCD(plot,'realJ1_ntrk12_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_width_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## ## plot = scaleQCD(plot,'realJ1_width12_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_MV1_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_jvf_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_ntrk_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## ## plot = scaleQCD(plot,'realJ2_ntrk12_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_width_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## ## plot = scaleQCD(plot,'realJ2_width12_Low',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_Low',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_sb_Low',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_lsb_Low',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'mllqq_hsb_Low',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
##############
plot = scaleQCD(plot,'chisquare_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'met',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
plot = scaleQCD(plot,'mll_High',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
plot = scaleQCD(plot,'ptll_High',chan,'incl',nmc[f_index]*tff_data.Integral()/tff_totmc.Integral(),nqcd[f_index]*tff_data.Integral()/tff_qcd.Integral())
plot = scaleQCD(plot,'mll_2j_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'leadleppt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'subleppt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'leadlepeta_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'sublepeta_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ1pt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ2pt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'corrJ1pt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'corrJ2pt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realZpt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'corrZpt_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realZm_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'corrZm_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'nJets_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'nbJets_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ1_MV1_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ1_jvf_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ1_ntrk_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_ntrk12_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ1_width_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ1_width12_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ2_MV1_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ2_jvf_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ2_ntrk_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_ntrk12_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'realJ2_width_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
## plot = scaleQCD(plot,'realJ2_width12_High',chan,'incl',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'mllqq_High',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'mllqq_sb_High',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'mllqq_lsb_High',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
plot = scaleQCD(plot,'mllqq_hsb_High',chan,'untag',nmc_2j[f_index]*tff_data_2j.Integral()/tff_totmc_2j.Integral(),nqcd_2j[f_index]*tff_data_2j.Integral()/tff_qcd_2j.Integral())
############################################################################################################



print '\nQCD extrapolation in mll SR:'
whichplot = 'mll_2j_High'
tagindex = 2 # inclusive
binmin = plot[QCD][whichplot][lepindex][tagindex].FindBin(70000)
binmax = plot[QCD][whichplot][lepindex][tagindex].FindBin(150000)
qcd_extr = plot[QCD][whichplot][lepindex][tagindex].Integral(binmin,binmax)
print 'QCD in 20-200 GeV is: %f \n\n' % (qcd_extr)

### scale to sidebands - test - ZSB
sb_incl_SFz = [0.,0.,0.]
sb_incl_SFz[lepindex] = calculateSFz(plot,'mllqq_sb_High',chan,'untag')

sb_low_SFz = [0.,0.,0.]
sb_low_SFz[lepindex] = calculateSFz(plot,'mllqq_lsb_High',chan,'untag')

sb_high_SFz = [0.,0.,0.]
sb_high_SFz[lepindex] = calculateSFz(plot,'mllqq_hsb_High',chan,'untag')

print '--------- Z+jets SF from the sidebands - ZSB -----------\n'
print ' inclusive SF = %f ' % sb_incl_SFz[lepindex]
print ' low mass SF = %f ' % sb_low_SFz[lepindex]
print ' high mass SF = %f ' % sb_high_SFz[lepindex]
print '----------------------------------------------------------'



lepindex = lepindexlist[0]
tagindex = 0 ## untag
#plot[Zlf_Alp]['mllqq_sb'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
#plot[Zlf_Alp]['mllqq_lsb'][lepindex][tagindex].Scale(sb_low_SFz[lepindex])
#plot[Zlf_Alp]['mllqq_hsb'][lepindex][tagindex].Scale(sb_high_SFz[lepindex])
##### plot[Zlf_Alp]['mllqq_Low'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
plot[Zlf_Alp]['mllqq_High'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
## plot[Zlf_Alp]['mll_exc'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
## plot[Zlf_Alp]['ptll_exc'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
## plot[Zlf_Alp]['mjj_exc'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
## plot[Zlf_Alp]['ptjj_exc'][lepindex][tagindex].Scale(sb_incl_SFz[lepindex])
####



for extsample in extsamplelist:
  for p in plotlist:
    for lepindex in lepindexlist:
      for tagindex in tagindexlist:
        if extsample.isBKG and not extsample.name=='totbkg' : plot[TotalBkg][p][lepindex][tagindex].Add(plot[extsample][p][lepindex][tagindex])




## plotinorder(plot,'chisquare',chan,'incl',False)
## plotinorder(plot,'mll',chan,'incl',False)
## dumpyields(plot,'mll',chan,'incl')
## plotinorder(plot,'mll_2j',chan,'incl',False)
## dumpyields(plot,'mll_2j',chan,'incl')
## plotinorder(plot,'mllqq',chan,'tag',False)
## dumpyields(plot,'mllqq',chan,'tag')
## plotinorder(plot,'mllqq',chan,'untag',False)
## dumpyields(plot,'mllqq',chan,'untag')
## plotinorder(plot,'leadleppt',chan,'incl',False)
## plotinorder(plot,'subleppt',chan,'incl',False)
## plotinorder(plot,'mllqq_sb',chan,'tag',False)
## dumpyields(plot,'mllqq_sb',chan,'tag')
## plotinorder(plot,'mllqq_sb',chan,'untag',False)
## dumpyields(plot,'mllqq_sb',chan,'untag')
## plotinorder(plot,'mjj',chan,'incl',False)
## plotinorder(plot,'mllqq_lsb',chan,'untag',False)
## dumpyields(plot,'mllqq_lsb',chan,'untag')
## plotinorder(plot,'mllqq_hsb',chan,'untag',False)
## dumpyields(plot,'mllqq_hsb',chan,'untag')
################################
## plotinorder(plot,'chisquare_Low',chan,'incl',False)
## plotinorder(plot,'mll_Low',chan,'incl',False)
## plotinorder(plot,'ptll_Low',chan,'incl',False)
## plotinorder(plot,'mll_2j_Low',chan,'incl',False)
## plotinorder(plot,'leadleppt_Low',chan,'incl',False)
## plotinorder(plot,'subleppt_Low',chan,'incl',False)
## plotinorder(plot,'leadlepeta_Low',chan,'incl',False)
## plotinorder(plot,'sublepeta_Low',chan,'incl',False)
## plotinorder(plot,'realJ1pt_Low',chan,'incl',False)
## plotinorder(plot,'realJ2pt_Low',chan,'incl',False)
## plotinorder(plot,'corrJ1pt_Low',chan,'incl',False)
## plotinorder(plot,'corrJ2pt_Low',chan,'incl',False)
## plotinorder(plot,'realZpt_Low',chan,'incl',False)
## plotinorder(plot,'corrZpt_Low',chan,'incl',False)
## plotinorder(plot,'realZm_Low',chan,'incl',False)
## plotinorder(plot,'corrZm_Low',chan,'incl',False)
## plotinorder(plot,'nJets_Low',chan,'incl',False)
## plotinorder(plot,'nbJets_Low',chan,'incl',False)
## plotinorder(plot,'realJ1_MV1_Low',chan,'incl',False)
## plotinorder(plot,'realJ1_jvf_Low',chan,'incl',False)
## plotinorder(plot,'realJ1_ntrk_Low',chan,'incl',False)
##     ## plotinorder(plot,'realJ1_ntrk12_Low',chan,'incl',False)
## plotinorder(plot,'realJ1_width_Low',chan,'incl',False)
##     ## plotinorder(plot,'realJ1_width12_Low',chan,'incl',False)
## plotinorder(plot,'realJ2_MV1_Low',chan,'incl',False)
## plotinorder(plot,'realJ2_jvf_Low',chan,'incl',False)
## plotinorder(plot,'realJ2_ntrk_Low',chan,'incl',False)
##     ## plotinorder(plot,'realJ2_ntrk12_Low',chan,'incl',False)
## plotinorder(plot,'realJ2_width_Low',chan,'incl',False)
##     ## plotinorder(plot,'realJ2_width12_Low',chan,'incl',False)
## plotinorder(plot,'mllqq_Low',chan,'untag',False)
## plotinorder(plot,'mllqq_sb_Low',chan,'untag',False)
## plotinorder(plot,'mllqq_lsb_Low',chan,'untag',False)
## plotinorder(plot,'mllqq_hsb_Low',chan,'untag',False)
##     ##############
plotinorder(plot,'chisquare_High',chan,'incl',False)
dumpyields(plot,'chisquare_High',chan,'incl')

plotinorder(plot,'met',chan,'incl',False)
dumpyields(plot,'met',chan,'incl')

plotinorder(plot,'mll_High',chan,'incl',False)
dumpyields(plot,'mll_High',chan,'incl')

plotinorder(plot,'ptll_High',chan,'incl',False)
dumpyields(plot,'ptll_High',chan,'incl')

plotinorder(plot,'mll_2j_High',chan,'incl',False)
dumpyields(plot,'mll_2j_High',chan,'incl')

plotinorder(plot,'leadleppt_High',chan,'incl',False)
dumpyields(plot,'leadleppt_High',chan,'incl')

plotinorder(plot,'subleppt_High',chan,'incl',False)
dumpyields(plot,'subleppt_High',chan,'incl')

plotinorder(plot,'leadlepeta_High',chan,'incl',False)
dumpyields(plot,'leadlepeta_High',chan,'incl')

plotinorder(plot,'sublepeta_High',chan,'incl',False)
dumpyields(plot,'sublepeta_High',chan,'incl')

plotinorder(plot,'realJ1pt_High',chan,'incl',False)
dumpyields(plot,'realJ1pt_High',chan,'incl')

plotinorder(plot,'realJ2pt_High',chan,'incl',False)
dumpyields(plot,'realJ2pt_High',chan,'incl')

plotinorder(plot,'corrJ1pt_High',chan,'incl',False)
dumpyields(plot,'corrJ1pt_High',chan,'incl')

plotinorder(plot,'corrJ2pt_High',chan,'incl',False)
dumpyields(plot,'corrJ2pt_High',chan,'incl')

plotinorder(plot,'realZpt_High',chan,'incl',False)
dumpyields(plot,'realZpt_High',chan,'incl')

plotinorder(plot,'corrZpt_High',chan,'incl',False)
dumpyields(plot,'corrZpt_High',chan,'incl')

plotinorder(plot,'realZm_High',chan,'incl',False)
dumpyields(plot,'realZm_High',chan,'incl')

plotinorder(plot,'corrZm_High',chan,'incl',False)
dumpyields(plot,'corrZm_High',chan,'incl')

plotinorder(plot,'nJets_High',chan,'incl',False)
dumpyields(plot,'nJets_High',chan,'incl')

plotinorder(plot,'nbJets_High',chan,'incl',False)
dumpyields(plot,'nbJets_High',chan,'incl')

plotinorder(plot,'realJ1_MV1_High',chan,'incl',False)
dumpyields(plot,'realJ1_MV1_High',chan,'incl')

plotinorder(plot,'realJ1_jvf_High',chan,'incl',False)
dumpyields(plot,'realJ1_jvf_High',chan,'incl')

plotinorder(plot,'realJ1_ntrk_High',chan,'incl',False)
dumpyields(plot,'realJ1_ntrk_High',chan,'incl')

## plotinorder(plot,'realJ1_ntrk12_High',chan,'incl',False)
## dumpyields(plot,'realJ1_ntrk12_High',chan,'incl')

plotinorder(plot,'realJ1_width_High',chan,'incl',False)
dumpyields(plot,'realJ1_width_High',chan,'incl')

plotinorder(plot,'realJ2_MV1_High',chan,'incl',False)
dumpyields(plot,'realJ2_MV1_High',chan,'incl')

plotinorder(plot,'realJ2_jvf_High',chan,'incl',False)
dumpyields(plot,'realJ2_jvf_High',chan,'incl')

plotinorder(plot,'realJ2_ntrk_High',chan,'incl',False)
dumpyields(plot,'realJ2_ntrk_High',chan,'incl')

## plotinorder(plot,'realJ2_ntrk12_High',chan,'incl',False)
## dumpyields(plot,'realJ2_ntrk12_High',chan,'incl')

plotinorder(plot,'realJ2_width_High',chan,'incl',False)
dumpyields(plot,'realJ2_width_High',chan,'incl')

## plotinorder(plot,'realJ2_width12_High',chan,'incl',False)
## dumpyields(plot,'realJ2_width12_High',chan,'incl')

plotinorder(plot,'mllqq_High',chan,'untag',False)
dumpyields(plot,'mllqq_High',chan,'incl')

plotinorder(plot,'mllqq_sb_High',chan,'untag',False)
dumpyields(plot,'mllqq_sb_High',chan,'incl')

plotinorder(plot,'mllqq_lsb_High',chan,'untag',False)
dumpyields(plot,'mllqq_lsb_High',chan,'incl')

plotinorder(plot,'mllqq_hsb_High',chan,'untag',False)
dumpyields(plot,'mllqq_hsb_High',chan,'incl')
##########################################################################################################



### save the histos for the limits
tagindex = 1 #tag
lepindex = lepindexlist[0]
fname = 'h2l2q_tag_'+chan+'_'+mass+'.root'
filetag = TFile(fname,'RECREATE')
for extsample in extsamplelist:
  if not extsample.name=='totbkg':
    plot[extsample]['mllqq_High'][lepindex][tagindex].SetName('mllqq_High_'+extsample.name)
    plot[extsample]['mllqq_High'][lepindex][tagindex].Write()
filetag.Write()
filetag.Close()

tagindex = 0 #untag
lepindex = lepindexlist[0]
fname = 'h2l2q_untag_'+chan+'_'+mass+'.root'
fileuntag = TFile(fname,'RECREATE')
for extsample in extsamplelist:
  if not extsample.name=='totbkg':
    plot[extsample]['mllqq_High'][lepindex][tagindex].SetName('mllqq_High_'+extsample.name)
    plot[extsample]['mllqq_High'][lepindex][tagindex].Write()
fileuntag.Write()
fileuntag.Close()
