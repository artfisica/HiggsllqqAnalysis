import glob

import sys
sys.path.append('./Higgs4lepAnalysis/python/resolution_studies')

def get_sample_filelist(path, name):
  import os

  os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep %s > __TMP_DPMLIST_RESOLUTIONSTUDIES__' % (path, name))

  result = []

  for line in open('__TMP_DPMLIST_RESOLUTIONSTUDIES__'):
    if name in line:
      result.append(line[:-1])

  os.system('rm  __TMP_DPMLIST_RESOLUTIONSTUDIES__')

  return result

# higgs class (container of mass, list of outputs in gg sample, list of outputs in vbf sample)
class higgsBoson:
  def __init__(self):
    self.mass = 0
    self.list_gg = []
    self.list_vbf = []
    self.min = 0;
    self.max = 0;
  def __init__(self, mass, list_gg, list_vbf, min, max):
    self.mass = mass
    self.list_gg = list_gg
    self.list_vbf = list_vbf
    if (min > max):
      self.max = min
      self.min = max
    else:
      self.min = min
      self.max = max

samples = [
#higgsBoson(110,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH110'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH110'),
#           80+(110-130),
#            150+(110-130)),
#higgsBoson(115,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH115'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH115'),
#           80+(115-130),
#            150+(115-130)),
#higgsBoson(120,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH120'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH120'),
#           80+(120-130),
#            150+(120-130)),
#higgsBoson(125,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH125'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH125'),
#           80+(125-130),
#            150+(125-130)),
#higgsBoson(130,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH130'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH130'),
#            80+(130-130),
#	     150+(130-130)),
 higgsBoson(130,
            glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*ggH130*/*root*'),
            [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH130'),
             80+(130-130),
 	     150+(130-130)),
#higgsBoson(135,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH135'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH135'),
#           80+(135-130),
#            150+(135-130)),
#higgsBoson(140,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH140'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH140'),
#           80+(140-130),
#            150+(140-130)),
#higgsBoson(145,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH145'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH145'),
#           80+(145-130),
#            150+(145-130)),
#higgsBoson(150,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH150'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH150'),
#           80+(150-130),
#            150+(150-130)),
#higgsBoson(180,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH180'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH180'),
#           130,
#            200),
#higgsBoson(260,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH260'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH260'),
#           170,
#            300),
#higgsBoson(360,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH360'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH360'),
#           250,
#            440),
#higgsBoson(460,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH460'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH460'),
#           250,
#            560),
#higgsBoson(600,
#           get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'ggH600'),
#           [],#get_sample_filelist('atlaslocalgroupdisk/user/vippolit/H4l000007/', 'VBFH600'),
#           350,
#            800),
]
