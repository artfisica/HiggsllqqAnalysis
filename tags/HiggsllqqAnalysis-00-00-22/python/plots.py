tagchan = {}
tagchan[0] = 'untag'
tagchan[1] = 'tag'
tagchan[2] = 'incl'

lepchan = {}
lepchan[0] = 'mumu'
lepchan[1] = 'mue'
lepchan[2] = 'ee'

plotlist = [
##     'chisquare_Low',
##     'mll_Low',
##     'ptll_Low',
##     'mll_2j_Low',
##     'leadleppt_Low',
##     'subleppt_Low',
##     'leadlepeta_Low',
##     'sublepeta_Low',
##     'realJ1pt_Low',
##     'realJ2pt_Low',
##     'corrJ1pt_Low',
##     'corrJ2pt_Low',
##     'realZpt_Low',
##     'corrZpt_Low',
##     'realZm_Low',
##     'corrZm_Low',
##     'nJets_Low',
##     'nbJets_Low',
##     'realJ1_MV1_Low',
##     'realJ1_jvf_Low',
##     'realJ1_ntrk_Low',
##     ## 'realJ1_ntrk12_Low',
##     'realJ1_width_Low',
##     ## 'realJ1_width12_Low',
##     'realJ2_MV1_Low',
##     'realJ2_jvf_Low',
##     'realJ2_ntrk_Low',
##     ## 'realJ2_ntrk12_Low',
##     'realJ2_width_Low',
##     ## 'realJ2_width12_Low',
##     'mllqq_Low',
##     'mllqq_sb_Low',
##     'mllqq_lsb_Low',
##     'mllqq_hsb_Low',
    ##############
    'chisquare_High',
    'met',
    'mll_High',
    'ptll_High',
    'mll_2j_High',
    'leadleppt_High',
    'subleppt_High',
    'leadlepeta_High',
    'sublepeta_High',
    'realJ1pt_High',
    'realJ2pt_High',
    'corrJ1pt_High',
    'corrJ2pt_High',
    'realZpt_High',
    'corrZpt_High',
    'realZm_High',
    'corrZm_High',
    'nJets_High',
    'nbJets_High',
    'realJ1_MV1_High',
    'realJ1_jvf_High',
    'realJ1_ntrk_High',
    ## 'realJ1_ntrk12_High',
    'realJ1_width_High',
    ## 'realJ1_width12_High',
    'realJ2_MV1_High',
    'realJ2_jvf_High',
    'realJ2_ntrk_High',
    ## 'realJ2_ntrk12_High',
    'realJ2_width_High',
    ## 'realJ2_width12_High',
    'mllqq_High',
    'mllqq_sb_High',
    'mllqq_lsb_High',
    'mllqq_hsb_High',
]

whattodraw = [
##     'chisquare',
##     'lepZ_m',
##     'lepZ_pt',
##     'lepZ_m',
##     'lep1_pt',
##     'lep2_pt',
##     'lep1_eta',
##     'lep2_eta',
##     'realJ1_pt',
##     'realJ2_pt',
##     'corrJ1_pt',
##     'corrJ2_pt',
##     'realZ_pt',
##     'corrZ_pt',
##     'realZ_m',
##     'corrZ_m',
##     'n_jets',
##     'n_b_jets',
##     'realJ1_MV1',
##     'realJ1_jvf',
##     'realJ1_ntrk',
##     'realJ1_ntrk12',
##     'realJ1_width',
##     'realJ1_width12',
##     'realJ2_MV1',
##     'realJ2_jvf',
##     'realJ2_ntrk',
##     'realJ2_ntrk12',
##     'realJ2_width',
##     'realJ2_width12',
##     'corrH_m',
##     'corrH_m',
##     'corrH_m',
##     'corrH_m',
    ##############
    'chisquare',
    'met',
    'lepZ_m',
    'lepZ_pt',
    'lepZ_m',
    'lep1_pt',
    'lep2_pt',
    'lep1_eta',
    'lep2_eta',
    'realJ1_pt',
    'realJ2_pt',
    'corrJ1_pt',
    'corrJ2_pt',
    'realZ_pt',
    'corrZ_pt',
    'realZ_m',
    'corrZ_m',
    'n_jets',
    'n_b_jets',
    'realJ1_MV1',
    'realJ1_jvf',
    'realJ1_ntrk',
    ## 'realJ1_ntrk12',
    'realJ1_width',
    ## 'realJ1_width12',
    'realJ2_MV1',
    'realJ2_jvf',
    'realJ2_ntrk',
    ## 'realJ2_ntrk12',
    'realJ2_width',
    ## 'realJ2_width12',
    'corrH_m',
    'corrH_m',
    'corrH_m',
    'corrH_m',
]




basecondition = [
##     ## 'chisquare_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && chisquare>=0',
##     ## 'mll_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0',
##     ## 'ptll_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0',
##     ## 'mll_2j_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && n_jets>=2',
##     ## 'leadleppt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0',
##     ## 'subleppt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000  && weight>=0',
##     ## 'leadlepeta_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0',
##     ## 'sublepeta_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0',
##     ## 'realJ1pt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ2pt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'corrJ1pt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && corrJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'corrJ2pt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && corrJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realZpt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && lepZ_m>20000 && lepZ_m<76000 && met<40000 && realJ1_pt>0 && realJ2_pt>0',
##     ## 'corrZpt_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && lepZ_m>20000 && lepZ_m<76000 && met<40000 && corrJ1_pt>0 && corrJ2_pt>0',
##     ## 'realZm_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && lepZ_m>20000 && lepZ_m<76000 && met<40000 && realJ1_pt>0 && realJ2_pt>0',
##     ## 'corrZm_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && lepZ_m>20000 && lepZ_m<76000 && met<40000 && corrJ1_pt>0 && corrJ2_pt>0',
##     ## 'nJets_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000',
##     ## 'nbJets_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000',
##     ## 'realJ1_MV1_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ1_jvf_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ1_ntrk_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ1_ntrk12_Low'
##     ## '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0',
##     ## 'realJ1_width_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ1_width12_Low'
##     ## '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ1_pt>0',
##     ## 'realJ2_MV1_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ2_jvf_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ2_ntrk_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ2_ntrk12_Low'
##     ## '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0',
##     ## 'realJ2_width_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0 && lepZ_m>20000 && lepZ_m<76000',
##     ## 'realJ2_width12_Low'
##     ## '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>15000 && weight>=0 && met<40000 && realJ2_pt>0',
##     ## 'mllqq_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && lepZ_m<76000 && chisquare>=0 && realZ_m>60000 && realZ_m<115000 && met<40000 && weight>=0', 
##     ## 'mllqq_sb_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && lepZ_m<76000 && chisquare>=0 && ((realZ_m>40000 && realZ_m<60000) || (realZ_m>115000 && realZ_m<160000)) && met<40000 && weight>=0',
##     ## 'mllqq_lsb_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && lepZ_m<76000 && chisquare>=0 && realZ_m>40000  && realZ_m<60000  && met<40000 && weight>=0',
##     ## 'mllqq_hsb_Low'
##     '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && lepZ_m<76000 && chisquare>=0 && realZ_m>115000 && realZ_m<160000 && met<40000 && weight>=0',
##     ###############################################################################################################
    ## 'chisquare_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0 && chisquare>=0',
    ## met
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'mll_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'ptll_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'mll_2j_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0 && n_jets>=2',
    ## 'leadleppt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'subleppt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000  && weight>=0',
    ## 'leadlepeta_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'sublepeta_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0',
    ## 'realJ1pt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ2pt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'corrJ1pt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && corrJ1_pt>0',
    ## 'corrJ2pt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && corrJ2_pt>0',
    ## 'realZpt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && lepZ_m>83000 && lepZ_m<99000 && met<60000 && realJ1_pt>0 && realJ2_pt>0',
    ## 'corrZpt_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && lepZ_m>83000 && lepZ_m<99000 && met<60000 && corrJ1_pt>0 && corrJ2_pt>0',
    ## 'realZm_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && lepZ_m>83000 && lepZ_m<99000 && met<60000 && realJ1_pt>0 && realJ2_pt>0',
    ## 'corrZm_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && lepZ_m>83000 && lepZ_m<99000 && met<60000 && corrJ1_pt>0 && corrJ2_pt>0',
    ## 'nJets_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0 && met<60000',
    ## 'nbJets_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>20000 && weight>=0 && met<60000',
    ## 'realJ1_MV1_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ1_jvf_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ1_ntrk_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ1_ntrk12_High'
    ## 'lepZ_m>20000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ1_width_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ1_width12_High'
    ## 'lepZ_m>20000 && weight>=0 && met<60000 && realJ1_pt>0',
    ## 'realJ2_MV1_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'realJ2_jvf_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'realJ2_ntrk_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'realJ2_ntrk12_High'
    ## 'lepZ_m>20000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'realJ2_width_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'realJ2_width12_High'
    ## 'lepZ_m>20000 && weight>=0 && met<60000 && realJ2_pt>0',
    ## 'mllqq_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && chisquare>=0 && realZ_m>70000 && realZ_m<105000 && met<60000 && weight>=0', 
    ## 'mllqq_sb_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && chisquare>=0 && ((realZ_m>40000 && realZ_m<70000) || (realZ_m>105000 && realZ_m<160000)) && met<60000 && weight>=0',
    ## 'mllqq_lsb_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && chisquare>=0 && realZ_m>40000  && realZ_m<70000  && met<60000 && weight>=0',
    ## 'mllqq_hsb_High'
    '((lep1_pt>20000 && lep2_pt>7000 && (trig_flag==1 || trig_flag==3)) || (lep1_pt>12000 && lep2_pt>12000 && channel==0 && (trig_flag==2 || trig_flag==3)) || (lep1_pt>14000 && lep2_pt>14000 && channel==2 && (trig_flag==2 || trig_flag==3)) || channel==1) && lepZ_m>83000 && lepZ_m<99000 && chisquare>=0 && realZ_m>105000 && realZ_m<160000 && met<60000 && weight>=0',
]




btagsf = [
##     ## 'chisquare_Low'
##     '',
##     ## 'mll_Low'
##     '',
##     ## 'ptll_Low'
##     '',
##     ## 'mll_2j_Low'
##     '*btagSF',
##     ## 'leadleppt_Low'
##     '',
##     ## 'subleppt_Low'
##     '',
##     ## 'leadlepeta_Low'
##     '',
##     ## 'sublepeta_Low'
##     '',
##     ## 'realJ1pt_Low'
##     '*btagSF',
##     ## 'realJ2pt_Low'
##     '*btagSF',
##     ## 'corrJ1pt_Low'
##     '*btagSF',
##     ## 'corrJ2pt_Low'
##     '*btagSF',
##     ## 'realZpt_Low'
##     '*btagSF',
##     ## 'corrZpt_Low'
##     '*btagSF',
##     ## 'realZm_Low'
##     '*btagSF',
##     ## 'corrZm_Low'
##     '*btagSF',
##     ## 'nJets_Low'
##     '*btagSF',
##     ## 'nbJets_Low'
##     '*btagSF',
##     ## 'realJ1_MV1_Low'
##     '*btagSF',
##     ## 'realJ1_jvf_Low'
##     '*btagSF',
##     ## 'realJ1_ntrk_Low'
##     '*btagSF',
##     ## 'realJ1_ntrk12_Low'
##     ## '*btagSF',
##     ## 'realJ1_width_Low'
##     '*btagSF',
##     ## 'realJ1_width12_Low'
##     ## '*btagSF',
##     ## 'realJ2_MV1_Low'
##     '*btagSF',
##     ## 'realJ2_jvf_Low'
##     '*btagSF',
##     ## 'realJ2_ntrk_Low'
##     '*btagSF',
##     ## 'realJ2_ntrk12_Low'
##     ## '*btagSF',
##     ## 'realJ2_width_Low'
##     '*btagSF',
##     ## 'realJ2_width12_Low'
##     ## '*btagSF',
##     ## 'mllqq_Low'
##     '*btagSF',
##     ## 'mllqq_sb_Low'
##     '*btagSF',
##     ## 'mllqq_lsb_Low'
##     '*btagSF',
##     ## 'mllqq_hsb_Low'
##     '*btagSF',
##     ##############
    ## 'chisquare_High'
    '',
    ## 'met'
    '',
    ## 'mll_High'
    '',
    ## 'ptll_High'
    '',
    ## 'mll_2j_High'
    '*btagSF',
    ## 'leadleppt_High'
    '',
    ## 'subleppt_High'
    '',
    ## 'leadlepeta_High'
    '',
    ## 'sublepeta_High'
    '',
    ## 'realJ1pt_High'
    '*btagSF',
    ## 'realJ2pt_High'
    '*btagSF',
    ## 'corrJ1pt_High'
    '*btagSF',
    ## 'corrJ2pt_High'
    '*btagSF',
    ## 'realZpt_High'
    '*btagSF',
    ## 'corrZpt_High'
    '*btagSF',
    ## 'realZm_High'
    '*btagSF',
    ## 'corrZm_High'
    '*btagSF',
    ## 'nJets_High'
    '*btagSF',
    ## 'nbJets_High'
    '*btagSF',
    ## 'realJ1_MV1_High'
    '*btagSF',
    ## 'realJ1_jvf_High'
    '*btagSF',
    ## 'realJ1_ntrk_High'
    '*btagSF',
    ## 'realJ1_ntrk12_High'
    ## '*btagSF',
    ## 'realJ1_width_High'
    '*btagSF',
    ## 'realJ1_width12_High'
    ## '*btagSF',
    ## 'realJ2_MV1_High'
    '*btagSF',
    ## 'realJ2_jvf_High'
    '*btagSF',
    ## 'realJ2_ntrk_High'
    '*btagSF',
    ## 'realJ2_ntrk12_High'
    ## '*btagSF',
    ## 'realJ2_width_High'
    '*btagSF',
    ## 'realJ2_width12_High'
    ## '*btagSF',
    ## 'mllqq_High'
    '*btagSF',
    ## 'mllqq_sb_High'
    '*btagSF',
    ## 'mllqq_lsb_High'
    '*btagSF',
    ## 'mllqq_hsb_High'
    '*btagSF',
]



nbins = [
##     ## 'chisquare_Low'
##     120,
##     ## 'mll_Low'
##     100,
##     ## 'ptll_Low'
##     100,
##     ## 'mll_2j_Low'
##     100,
##     ## 'leadleppt_Low'
##     100,
##     ## 'subleppt_Low'
##     100,
##     ## 'leadlepeta_Low'
##     10,
##     ## 'sublepeta_Low'
##     10,
##     ## 'realJ1pt_Low'
##     100,
##     ## 'realJ2pt_Low'
##     100,
##     ## 'corrJ1pt_Low'
##     100,
##     ## 'corrJ2pt_Low'
##     100,
##     ## 'realZpt_Low'
##     100,
##     ## 'corrZpt_Low'
##     100,
##     ## 'realZm_Low'
##     100,
##     ## 'corrZm_Low'
##     100,
##     ## 'nJets_Low'
##     10,
##     ## 'nbJets_Low'
##     10,
##     ## 'realJ1_MV1_Low'
##     100,
##     ## 'realJ1_jvf_Low'
##     100,
##     ## 'realJ1_ntrk_Low'
##     20,
##     ## 'realJ1_ntrk12_Low'
##     ## 20,
##     ## 'realJ1_width_Low'
##     100,
##     ## 'realJ1_width12_Low'
##     ## 100,
##     ## 'realJ2_MV1_Low'
##     100,
##     ## 'realJ2_jvf_Low'
##     100,
##     ## 'realJ2_ntrk_Low'
##     20,
##     ## 'realJ2_ntrk12_Low'
##     ## 20,
##     ## 'realJ2_width_Low'
##     100,
##     ## 'realJ2_width12_Low'
##     ## 100,
##     ## 'mllqq_Low'
##     200,
##     ## 'mllqq_sb_Low'
##     200,
##     ## 'mllqq_lsb_Low'
##     200,
##     ## 'mllqq_hsb_Low'
##     200,
##     ##############
    ## 'chisquare_High'
    120,
    ## 'met'
    150,
    ## 'mll_High'
    100,
    ## 'ptll_High'
    100,
    ## 'mll_2j_High'
    100,
    ## 'leadleppt_High'
    100,
    ## 'subleppt_High'
    100,
    ## 'leadlepeta_High'
    10,
    ## 'sublepeta_High'
    10,
    ## 'realJ1pt_High'
    100,
    ## 'realJ2pt_High'
    100,
    ## 'corrJ1pt_High'
    100,
    ## 'corrJ2pt_High'
    100,
    ## 'realZpt_High'
    100,
    ## 'corrZpt_High'
    100,
    ## 'realZm_High'
    100,
    ## 'corrZm_High'
    100,
    ## 'nJets_High'
    10,
    ## 'nbJets_High'
    10,
    ## 'realJ1_MV1_High'
    100,
    ## 'realJ1_jvf_High'
    100,
    ## 'realJ1_ntrk_High'
    20,
    ## 'realJ1_ntrk12_High'
    ## 20,
    ## 'realJ1_width_High'
    100,
    ## 'realJ1_width12_High'
    ## 100,
    ## 'realJ2_MV1_High'
    100,
    ## 'realJ2_jvf_High'
    100,
    ## 'realJ2_ntrk_High'
    20,
    ## 'realJ2_ntrk12_High'
    ## 20,
    ## 'realJ2_width_High'
    100,
    ## 'realJ2_width12_High'
    ## 100,
    ## 'mllqq_High'
    200,
    ## 'mllqq_sb_High'
    200,
    ## 'mllqq_lsb_High'
    200,
    ## 'mllqq_hsb_High'
    200,
]



binmin = [
##     ## 'chisquare_Low'
##     0,
##     ## 'mll_Low'
##     0,
##     ## 'ptll_Low'
##     0,
##     ## 'mll_2j_Low'
##     0,
##     ## 'leadleppt_Low'
##     0,
##     ## 'subleppt_Low'
##     0,
##     ## 'leadlepeta_Low'
##     0,
##     ## 'sublepeta_Low'
##     0,
##     ## 'realJ1pt_Low'
##     0,
##     ## 'realJ2pt_Low'
##     0,
##     ## 'corrJ1pt_Low'
##     0,
##     ## 'corrJ2pt_Low'
##     0,
##     ## 'realZpt_Low'
##     0,
##     ## 'corrZpt_Low'
##     0,
##     ## 'realZm_Low'
##     10000,
##     ## 'corrZm_Low'
##     10000,
##     ## 'nJets_Low'
##     0,
##     ## 'nbJets_Low'
##     0,
##     ## 'realJ1_MV1_Low'
##     -1,
##     ## 'realJ1_jvf_Low'
##     -1,
##     ## 'realJ1_ntrk_Low'
##     0,
##     ## 'realJ1_ntrk12_Low'
##     ## 0,
##     ## 'realJ1_width_Low'
##     0,
##     ## 'realJ1_width12_Low'
##     ## 0,
##     ## 'realJ2_MV1_Low'
##     -1,
##     ## 'realJ2_jvf_Low'
##     -1,
##     ## 'realJ2_ntrk_Low'
##     0,
##     ## 'realJ2_ntrk12_Low'
##     ## 0,
##     ## 'realJ2_width_Low'
##     0,
##     ## 'realJ2_width12_Low'
##     ## 0,
##     ## 'mllqq_Low'
##     20000,
##     ## 'mllqq_sb_Low'
##     20000,
##     ## 'mllqq_lsb_Low'
##     20000,
##     ## 'mllqq_hsb_Low'
##     20000,
##     ##############
    ## 'chisquare_High'
    0,
    ## 'met'
    0,
    ## 'mll_High'
    0,
    ## 'ptll_High'
    0,
    ## 'mll_2j_High'
    0,
    ## 'leadleppt_High'
    0,
    ## 'subleppt_High'
    0,
    ## 'leadlepeta_High'
    0,
    ## 'sublepeta_High'
    0,
    ## 'realJ1pt_High'
    0,
    ## 'realJ2pt_High'
    0,
    ## 'corrJ1pt_High'
    0,
    ## 'corrJ2pt_High'
    0,
    ## 'realZpt_High'
    0,
    ## 'corrZpt_High'
    0,
    ## 'realZm_High'
    10000,
    ## 'corrZm_High'
    10000,
    ## 'nJets_High'
    0,
    ## 'nbJets_High'
    0,
    ## 'realJ1_MV1_High'
    -1,
    ## 'realJ1_jvf_High'
    -1,
    ## 'realJ1_ntrk_High'
    0,
    ## 'realJ1_ntrk12_High'
    ## 0,
    ## 'realJ1_width_High'
    0,
    ## 'realJ1_width12_High'
    ## 0,
    ## 'realJ2_MV1_High'
    -1,
    ## 'realJ2_jvf_High'
    -1,
    ## 'realJ2_ntrk_High'
    0,
    ## 'realJ2_ntrk12_High'
    ## 0,
    ## 'realJ2_width_High'
    0,
    ## 'realJ2_width12_High'
    ## 0,
    ## 'mllqq_High'
    20000,
    ## 'mllqq_sb_High'
    20000,
    ## 'mllqq_lsb_High'
    20000,
    ## 'mllqq_hsb_High'
    20000,
]




binmax = [
##     ## 'chisquare_Low'
##     30,
##     ## 'mll_Low'
##     200000,
##     ## 'ptll_Low'
##     200000,
##     ## 'mll_2j_Low'
##     200000,
##     ## 'leadleppt_Low'
##     150000,
##     ## 'subleppt_Low'
##     150000,
##     ## 'leadlepeta_Low'
##     2.5,
##     ## 'sublepeta_Low'
##     2.5,
##     ## 'realJ1pt_Low'
##     200000,
##     ## 'realJ2pt_Low'
##     200000,
##     ## 'corrJ1pt_Low'
##     200000,
##     ## 'corrJ2pt_Low'
##     200000,
##     ## 'realZpt_Low'
##     200000,
##     ## 'corrZpt_Low'
##     200000,
##     ## 'realZm_Low'
##     200000,
##     ## 'corrZm_Low'
##     200000,
##     ## 'nJets_Low'
##     10,
##     ## 'nbJets_Low'
##     10,
##     ## 'realJ1_MV1_Low'
##     1,
##     ## 'realJ1_jvf_Low'
##     1,
##     ## 'realJ1_ntrk_Low'
##     20,
##     ## 'realJ1_ntrk12_Low'
##     ## 20,
##     ## 'realJ1_width_Low'
##     1,
##     ## 'realJ1_width12_Low'
##     ## 1,
##     ## 'realJ2_MV1_Low'
##     1,
##     ## 'realJ2_jvf_Low'
##     1,
##     ## 'realJ2_ntrk_Low'
##     20,
##     ## 'realJ2_ntrk12_Low'
##     ## 20,
##     ## 'realJ2_width_Low'
##     1,
##     ## 'realJ2_width12_Low'
##     ## 1,
##     ## 'mllqq_Low'
##     1000000,
##     ## 'mllqq_sb_Low'
##     1000000,
##     ## 'mllqq_lsb_Low'
##     1000000,
##     ## 'mllqq_hsb_Low'
##     1000000,
##     ##############
    ## 'chisquare_High'
    30,
    ## 'met'
    150000,
    ## 'mll_High'
    200000,
    ## 'ptll_High'
    200000,
    ## 'mll_2j_High'
    200000,
    ## 'leadleppt_High'
    150000,
    ## 'subleppt_High'
    150000,
    ## 'leadlepeta_High'
    2.5,
    ## 'sublepeta_High'
    2.5,
    ## 'realJ1pt_High'
    200000,
    ## 'realJ2pt_High'
    200000,
    ## 'corrJ1pt_High'
    200000,
    ## 'corrJ2pt_High'
    200000,
    ## 'realZpt_High'
    200000,
    ## 'corrZpt_High'
    200000,
    ## 'realZm_High'
    200000,
    ## 'corrZm_High'
    200000,
    ## 'nJets_High'
    10,
    ## 'nbJets_High'
    10,
    ## 'realJ1_MV1_High'
    1,
    ## 'realJ1_jvf_High'
    1,
    ## 'realJ1_ntrk_High'
    20,
    ## 'realJ1_ntrk12_High'
    ## 20,
    ## 'realJ1_width_High'
    1,
    ## 'realJ1_width12_High'
    ## 1,
    ## 'realJ2_MV1_High'
    1,
    ## 'realJ2_jvf_High'
    1,
    ## 'realJ2_ntrk_High'
    20,
    ## 'realJ2_ntrk12_High'
    ## 20,
    ## 'realJ2_width_High'
    1,
    ## 'realJ2_width12_High'
    ## 1,
    ## 'mllqq_High'
    1000000,
    ## 'mllqq_sb_High'
    1000000,
    ## 'mllqq_lsb_High'
    1000000,
    ## 'mllqq_hsb_High'
    1000000,
]
