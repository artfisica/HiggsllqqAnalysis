### dumps full details for a list of candidates
### author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *

import os, re
import glob

import sys
sys.path.append('./datasets')

interesting_candidates = [
[182284,   91584073],
[182486,   33852510],
[182766,    5404925],
[183003,   44433120],
[183003,  121099951],
[183045,   28329677],
[183081,  121479214],
[183391,   19834577],
[183426,   47756740],
[183602,     282919],
[184130,  194694606],
[186156,   65491657],
[186361,   29706533],
[187219,   88203394],
[187811,   40399631],
[189207,   79774710],
[189280,   82801561],
[189280,  128083498],
[189280,  143576946],
[189561,   20659041],
[189561,   48300215],
[189693,   10714212],
[189822,   75634934],
[190116,   60445481],
[190116,   63457816],
[190256,  175308437],
[190643,   23010890],
[190872,   52781235],
[191190,    9522152],
[191426,   60906769],
[191513,    9086448],
[191635,   48018112],
[191676,    1888359],
[182747,  136528490],
[182787,   35518831],
[183216,   75692579],
[184022,   20046902],
[184022,   78541915],
[186216,   10253640],
[186721,   22103472],
[186729,  203362752],
[186934,   65787798],
[187219,   24256525],
[187219,   34502707],
[187453,   34960141],
[187552,    3744932],
[187763,   83732606],
[189483,   33468656],
[189751,   51800361],
[190256,   56600537],
[190933,   99272087],
[191138,   17388332],
[191139,    5871977],
[191635,    2200900],
[191920,    7775570],
[179710,   25946709],
[180636,   71391739],
[180710,   37143864],
[182486,   21528951],
[182747,   63217197],
[182796,   74566644],
[183407,  136901836],
[183462,   75344317],
[186877,   12509901],
[186923,   96974859],
[189483,    1021987],
[189719,   37988693],
[190046,    8638208],
[190300,   17344710],
[190300,  121067450],
[190878,   50034828],
[190975,   20471852],
[191138,   15762515],
[191150,    5742674],
[191150,   45707611],
[191190,   76273161],
[183426,   50303812],
[186399,   14250520],
[186877,   84622334],
[187014,  105211056],
[189242,    7233912],
[189561,  105481981],
[189781,    8619753],
[190878,   14892058],
[190878,   57044890],
[190975,   62905396],
[191218,    1072214],
[191428,   25718643],
[200926,   17712380],
[201494,   26817712],
[202668,   26299894],
[203195,    5763465],
[203258,   59505222],
[203277,   22949094],
[203353,   26647061],
[203432,   46880033],
[203456,   24914034],
[203524,   62321499],
[203602,   80248134],
[203636,   34504484],
[203745,   63994592],
[203779,   46710128],
[203876,    9282788],
[203934,  100600041],
[204026,   87732822],
[204071,   30963264],
[204073,   18405446],
[204153,   22992989],
[204153,   32436165],
[204158,  108769641],
[204240,   52902065],
[204416,   20491680],
[204564,   25416035],
[204564,  149682166],
[204763,   64671324],
[204763,  198344978],
[204769,   71902630],
[204769,   82599793],
[204769,   84802829],
[204976,   52368897],
[205016,   37425673],
[205017,    7255669],
[205071,   55053974],
[205071,  222074238],
[205112,   29305779],
[205112,   58173349],
[201113,   36934150],
[202991,   29492591],
[203335,   44927851],
[203432,   58320601],
[203524,   64340754],
[203602,   82614360],
[203602,   99681996],
[203602,  106701658],
[203636,   95484441],
[203680,   44088690],
[203719,   20441551],
[203739,   73042609],
[203934,   89714329],
[203934,   91585470],
[204240,   29054747],
[204265,   11531490],
[204265,   15870479],
[204564,   72798274],
[204564,   76786292],
[204564,  193961472],
[204564,  228171673],
[204633,    4809722],
[204668,   60678687],
[204726,   20583917],
[204910,   22993546],
[204932,   53150539],
[204932,   60690285],
[204954,   19454464],
[204955,   90366740],
[205055,   13546817],
[205071,  193422310],
[205112,   46981864],
[200987,   50726675],
[201113,   29045106],
[201113,   34945963],
[201138,   27887257],
[201190,   13200815],
[202798,   43736485],
[203027,   71454016],
[203195,   25217284],
[203228,   15214720],
[203258,  105740575],
[203353,   79978481],
[203432,   20659622],
[203523,   13350735],
[203636,   10916804],
[203636,   69879354],
[203739,   67630085],
[203876,   96963848],
[204153,   33235991],
[204158,   18242755],
[204240,   31409138],
[204240,  151358568],
[204265,   72546153],
[204416,   42212637],
[204474,   56451245],
[204474,   89082513],
[204763,    5630931],
[204763,   67344918],
[204910,   95498501],
[204955,   16392574],
[204955,   72469686],
[205112,   25095415],
[205113,   12611816],
[201113,   55604596],
[201257,   21308431],
[201257,  121862343],
[201289,   10043905],
[201556,    6395332],
[201556,   10787809],
[202798,    7203716],
[203258,  114413312],
[203335,   49048571],
[203353,    7225045],
[203353,   28661403],
[203524,   28303221],
[203636,   65210779],
[203680,   10638148],
[203719,   12262192],
[203745,   75039428],
[203745,   90860257],
[203876,   36135200],
[203876,   60324400],
[203934,   70728810],
[203934,   98492930],
[204071,   32311568],
[204073,   38389153],
[204153,   26755807],
[204265,  223273173],
[204474,  131005091],
[204474,  134814016],
[204564,   68230882],
[204564,  198222712],
[204668,  141229289],
[204763,   95056361],
[204769,   38987331],
[204772,   81936527],
[204932,   57770438],
[204954,   36015359],
[205017,   18590216],
[205055,  201414970],
[205071,   36046056],
[205112,    5259842],
[205113,   28375316],
#[204769, 82599793],
#[204769, 71902630],
# 2012, ~900 GeV
#[204564, 68230882],
# 2012, [120-130] GeV
#[204763,  95056361],
#[204153,  22992989],
#[205113,  12611816],
#[204769,  82599793],
#[204769,  71902630],
#[203602,  82614360],
#[204910,  22993546],
#[205071,  36046056],
#[204564,  25416035],
#[203258, 105740575],
]

#output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*physics*/*root*')
output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_201*_v20/*physics*/*root*')

for filename in output_list:
 #print 'reading ' + filename
  file = TFile(filename)
  tree = file.Get("candidates")

  for i in range(tree.GetEntries()):
    tree.GetEntry(i)

    for candidate in interesting_candidates:
      if (tree.run == candidate[0] and tree.event == candidate[1]):
	print ''
	print ' [from %s]' % (filename)
	print '##########################################################'
        print 'run %d event %d (lbn %d)' % (tree.run, tree.event, tree.lbn)
	print 'candidate of type %d passes cut %d and has selected == %d and closesttozz == %d' % (tree.type, tree.last, tree.selected, tree.closesttozz)
	print '  Higgs mass = %lf (constrained = %lf)' % (tree.H_m, tree.H_m_constrained)
	print '     Z1 mass = %lf' % (tree.Z1_m)
	print '        -> lepplus:'
	print '                    pt = %lf' % (tree.Z1_lepplus_pt)
	print '                   eta = %lf' % (tree.Z1_lepplus_eta)
	print '                   phi = %lf' % (tree.Z1_lepplus_phi)
	print '                     m = %lf' % (tree.Z1_lepplus_m)
	print '                author = %d'  % (tree.Z1_lepplus_author)
	print '                  isSA = %d'  % (tree.Z1_lepplus_isSA)
	print '              ptcone20 = %lf' % (tree.Z1_lepplus_ptcone20)
	print '              etcone20 = %lf' % (tree.Z1_lepplus_etcone20)
	print '                    d0 = %lf' % (tree.Z1_lepplus_d0)
	print '                d0_sig = %lf' % (tree.Z1_lepplus_d0_sig)
	print '     ptcone20_final/pt = %lf' % (tree.Z1_lepplus_ptcone20_final / tree.Z1_lepplus_pt)
	print '     etcone20_final/pt = %lf' % (tree.Z1_lepplus_etcone20_final / tree.Z1_lepplus_pt)
	print '             d0/d0_sig = %lf' % (tree.Z1_lepplus_d0/tree.Z1_lepplus_d0_sig)
	print '        -> lepminus:'
	print '                    pt = %lf' % (tree.Z1_lepminus_pt)
	print '                   eta = %lf' % (tree.Z1_lepminus_eta)
	print '                   phi = %lf' % (tree.Z1_lepminus_phi)
	print '                     m = %lf' % (tree.Z1_lepminus_m)
	print '                author = %d'  % (tree.Z1_lepminus_author)
	print '                  isSA = %d'  % (tree.Z1_lepminus_isSA)
	print '              ptcone20 = %lf' % (tree.Z1_lepminus_ptcone20)
	print '              etcone20 = %lf' % (tree.Z1_lepminus_etcone20)
	print '                    d0 = %lf' % (tree.Z1_lepminus_d0)
	print '                d0_sig = %lf' % (tree.Z1_lepminus_d0_sig)
	print '     ptcone20_final/pt = %lf' % (tree.Z1_lepminus_ptcone20_final / tree.Z1_lepminus_pt)
	print '     etcone20_final/pt = %lf' % (tree.Z1_lepminus_etcone20_final / tree.Z1_lepminus_pt)
	print '             d0/d0_sig = %lf' % (tree.Z1_lepminus_d0/tree.Z1_lepminus_d0_sig)
	print '     Z2 mass = %lf' % (tree.Z2_m)
	print '        -> lepplus:'
	print '                    pt = %lf' % (tree.Z2_lepplus_pt)
	print '                   eta = %lf' % (tree.Z2_lepplus_eta)
	print '                   phi = %lf' % (tree.Z2_lepplus_phi)
	print '                     m = %lf' % (tree.Z2_lepplus_m)
	print '                author = %d'  % (tree.Z2_lepplus_author)
	print '                  isSA = %d'  % (tree.Z2_lepplus_isSA)
	print '              ptcone20 = %lf' % (tree.Z2_lepplus_ptcone20)
	print '              etcone20 = %lf' % (tree.Z2_lepplus_etcone20)
	print '                    d0 = %lf' % (tree.Z2_lepplus_d0)
	print '                d0_sig = %lf' % (tree.Z2_lepplus_d0_sig)
	print '     ptcone20_final/pt = %lf' % (tree.Z2_lepplus_ptcone20_final / tree.Z2_lepplus_pt)
	print '     etcone20_final/pt = %lf' % (tree.Z2_lepplus_etcone20_final / tree.Z2_lepplus_pt)
	print '             d0/d0_sig = %lf' % (tree.Z2_lepplus_d0/tree.Z2_lepplus_d0_sig)
	print '        -> lepminus:'
	print '                    pt = %lf' % (tree.Z2_lepminus_pt)
	print '                   eta = %lf' % (tree.Z2_lepminus_eta)
	print '                   phi = %lf' % (tree.Z2_lepminus_phi)
	print '                     m = %lf' % (tree.Z2_lepminus_m)
	print '                author = %d'  % (tree.Z2_lepminus_author)
	print '                  isSA = %d'  % (tree.Z2_lepminus_isSA)
	print '              ptcone20 = %lf' % (tree.Z2_lepminus_ptcone20)
	print '              etcone20 = %lf' % (tree.Z2_lepminus_etcone20)
	print '                    d0 = %lf' % (tree.Z2_lepminus_d0)
	print '                d0_sig = %lf' % (tree.Z2_lepminus_d0_sig)
	print '     ptcone20_final/pt = %lf' % (tree.Z2_lepminus_ptcone20_final / tree.Z2_lepminus_pt)
	print '     etcone20_final/pt = %lf' % (tree.Z2_lepminus_etcone20_final / tree.Z2_lepminus_pt)
	print '             d0/d0_sig = %lf' % (tree.Z2_lepminus_d0/tree.Z2_lepminus_d0_sig)

        Z1_lepplus_4m = TLorentzVector()
        Z1_lepminus_4m = TLorentzVector()
        Z2_lepplus_4m = TLorentzVector()
        Z2_lepminus_4m = TLorentzVector()

        Z1_lepplus_4m.SetPtEtaPhiM(tree.Z1_lepplus_pt, tree.Z1_lepplus_eta, tree.Z1_lepplus_phi, tree.Z1_lepplus_m)
        Z1_lepminus_4m.SetPtEtaPhiM(tree.Z1_lepminus_pt, tree.Z1_lepminus_eta, tree.Z1_lepminus_phi, tree.Z1_lepminus_m)
        Z2_lepplus_4m.SetPtEtaPhiM(tree.Z2_lepplus_pt, tree.Z2_lepplus_eta, tree.Z2_lepplus_phi, tree.Z2_lepplus_m)
        Z2_lepminus_4m.SetPtEtaPhiM(tree.Z2_lepminus_pt, tree.Z2_lepminus_eta, tree.Z2_lepminus_phi, tree.Z2_lepminus_m)

	quadrimomenta = [Z1_lepplus_4m, Z1_lepminus_4m, Z2_lepplus_4m, Z2_lepminus_4m]

	for i in range(len(quadrimomenta)):
	  for j in range(len(quadrimomenta)):
	    print '        deltaR [%d, %d] is %lf ' % (i, j, quadrimomenta[i].DeltaR(quadrimomenta[j]))
	
	print 'comes from ' + tree.filename[0]
