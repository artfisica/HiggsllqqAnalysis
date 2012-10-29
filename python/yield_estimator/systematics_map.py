# adapted from table 43 of https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG2/ConfNotes/NoteZZ4linternal_March2012/systematics.tex
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import uncertainty_computer as uc

syst_unc_mass = {}
#                        ELECTRON                   MUON 
#              mH   4mu  2e2mu  2mu2e   4e     4mu   2e2mu  2mu2e  4e 
syst_unc_mass[110] = [0, 2.796, 5.698, 7.111, 0.161, 0.112, 0.102, 0]
syst_unc_mass[115] = [0, 2.685, 6.770, 7.622, 0.162, 0.114, 0.103, 0]
syst_unc_mass[120] = [0, 2.710, 6.911, 7.898, 0.167, 0.152, 0.105, 0]
syst_unc_mass[125] = [0, 2.723, 7.040, 7.796, 0.187, 0.126, 0.108, 0]
syst_unc_mass[130] = [0, 2.739, 6.581, 7.399, 0.191, 0.151, 0.107, 0]
syst_unc_mass[140] = [0, 2.795, 5.777, 6.212, 0.175, 0.137, 0.123, 0]
syst_unc_mass[145] = [0, 2.829, 5.278, 5.819, 0.185, 0.130, 0.107, 0]
syst_unc_mass[150] = [0, 2.813, 4.838, 5.423, 0.195, 0.125, 0.110, 0]
syst_unc_mass[160] = [0, 2.830, 4.256, 4.977, 0.175, 0.119, 0.109, 0]
syst_unc_mass[170] = [0, 2.791, 3.686, 4.392, 0.170, 0.105, 0.112, 0]
syst_unc_mass[180] = [0, 2.718, 3.192, 3.902, 0.166, 0.117, 0.109, 0]
syst_unc_mass[185] = [0, 2.666, 3.080, 3.730, 0.163, 0.107, 0.106, 0]
syst_unc_mass[190] = [0, 2.691, 3.060, 3.729, 0.164, 0.113, 0.115, 0]
syst_unc_mass[195] = [0, 2.683, 3.039, 3.715, 0.168, 0.109, 0.112, 0]
syst_unc_mass[200] = [0, 2.690, 3.032, 3.718, 0.167, 0.113, 0.111, 0]
syst_unc_mass[210] = [0, 2.737, 3.109, 3.779, 0.176, 0.118, 0.110, 0]
syst_unc_mass[240] = [0, 2.743, 2.977, 3.708, 0.195, 0.120, 0.114, 0]
syst_unc_mass[260] = [0, 2.674, 2.973, 3.625, 0.183, 0.130, 0.121, 0]
syst_unc_mass[280] = [0, 2.641, 2.884, 3.550, 0.192, 0.118, 0.130, 0]
syst_unc_mass[300] = [0, 2.632, 2.886, 3.502, 0.186, 0.125, 0.124, 0]
syst_unc_mass[320] = [0, 2.593, 2.819, 3.447, 0.188, 0.123, 0.116, 0]
syst_unc_mass[340] = [0, 2.571, 2.742, 3.380, 0.194, 0.122, 0.126, 0]
syst_unc_mass[360] = [0, 2.563, 2.728, 3.293, 0.189, 0.126, 0.122, 0]
syst_unc_mass[380] = [0, 2.549, 2.701, 3.304, 0.189, 0.129, 0.126, 0]
syst_unc_mass[420] = [0, 2.530, 2.672, 3.220, 0.198, 0.127, 0.126, 0]
syst_unc_mass[440] = [0, 2.494, 2.636, 3.156, 0.194, 0.124, 0.122, 0]
syst_unc_mass[460] = [0, 2.512, 2.624, 3.175, 0.195, 0.132, 0.118, 0]
syst_unc_mass[480] = [0, 2.486, 2.597, 3.168, 0.190, 0.123, 0.128, 0]
syst_unc_mass[500] = [0, 2.491, 2.619, 3.113, 0.186, 0.124, 0.125, 0]
syst_unc_mass[520] = [0, 2.436, 2.618, 3.071, 0.189, 0.125, 0.124, 0]
syst_unc_mass[540] = [0, 2.455, 2.610, 3.077, 0.186, 0.131, 0.125, 0]
syst_unc_mass[560] = [0, 2.449, 2.595, 3.015, 0.193, 0.130, 0.125, 0]
syst_unc_mass[580] = [0, 2.447, 2.573, 3.035, 0.193, 0.131, 0.125, 0]
syst_unc_mass[600] = [0, 2.417, 2.585, 3.000, 0.190, 0.126, 0.119, 0]


def get_syst_unc(mass, chan):
  if (chan == '4mu'):
    return uc.sum_sq([syst_unc_mass[mass][0]/100., syst_unc_mass[mass][4]]/100.)
  elif (chan == '2e2mu'):
    return uc.sum_sq([syst_unc_mass[mass][1]/100., syst_unc_mass[mass][5]]/100.)
  elif (chan == '2mu2e'):
    return uc.sum_sq([syst_unc_mass[mass][2]/100., syst_unc_mass[mass][6]]/100.)
  elif (chan == '4e'):
    return uc.sum_sq([syst_unc_mass[mass][3]/100., syst_unc_mass[mass][7]]/100.)
  else:
    print 'ERROR: called get_syst_unc(%d, %s)' % (mass, chan)
    return 0
