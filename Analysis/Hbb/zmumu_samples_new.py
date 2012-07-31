import ROOT
SAMPLES={
     # key :           (         rootfile,                                   events, cross section)
     "DATA"           :("2012hsts/Hbb_Zmumu__SingleMuon2012_all.root",                0. ,     0.0),
     "HIGGS"          :("2012hsts/Hbb_Zmumu__HToBB_M-125.root"           ,            879458 ,   22.97),
     "ZJets_udscg_ptL":("2012hsts/Hbb_Zmumu__DYJetsToLL_pT70-100_udscgJets.root",1.41339e+06 , 62140.0),
     "ZJets_bJets_ptL":("2012hsts/Hbb_Zmumu__DYJetsToLL_pT70-100_bJets.root"    ,1.41339e+06 , 62140.0),
     "ZJets_udscg_ptH":("2012hsts/Hbb_Zmumu__DYJetsToLL_pT100_udscgJets.root"   ,2.64e+06    , 40510.0),
     "ZJets_bJets_ptH":("2012hsts/Hbb_Zmumu__DYJetsToLL_pT100_bJets.root"       ,2.64e+06    , 40510.0),
     "ZZ"             :("2012hsts/Hbb_Zmumu__ZZ.root"                           ,4.429895e+06, 8255.61),
     "WW"             :("2012hsts/Hbb_Zmumu__WW.root"                           ,4.341143e+06, 57109.7),
     "WZ"             :("2012hsts/Hbb_Zmumu__WZ.root"                           ,3.764351e+06, 32316.1),
     "TTbar"          :("2012hsts/Hbb_Zmumu__TTJets.root"                       ,6.27611e+06,  225197.),
     "ST-s"           :("2012hsts/Hbb_Zmumu__T_s-channel.root"                  ,2.39999e+05,  3893.94),
     "STbar-s"        :("2012hsts/Hbb_Zmumu__Tbar_s-channel.root"               ,1.39973e+05,  1757.76),
     ## "ST-t"           :("2012hsts/Hbb_Zmumu__T_t-channel.root"                  ,xx,           30004.2),
     "STbar-t"        :("2012hsts/Hbb_Zmumu__Tbar_t-channel.root"               ,1.455067e+06,  55531.0),
     "ST-tW"          :("2012hsts/Hbb_Zmumu__T_tW-channel-DR.root"              ,4.37656e+05,  11177.3),
     "STbar-tW"       :("2012hsts/Hbb_Zmumu__Tbar_tW-channel-DR.root"           ,4.63458e+05,  11177.3),
     }

##SAMPLES={
##    # key :           (         rootfile,                                   events, cross section)
##    "DATA"           :("cHists/Hbb_Zmumu__SingleMuon2012_all.root",                0. ,     0.0),
##    "HIGGS"          :("cHists/Hbb_Zmumu__HToBB_M-120.root"           ,            950000 ,   22.97),
##    "ZJets_udscg_ptL":("cHists/Hbb_Zmumu__DYJetsToLL_pT70-100_udscgJets.root",1.41339e+06 , 62140.0),
##    "ZJets_bJets_ptL":("cHists/Hbb_Zmumu__DYJetsToLL_pT70-100_bJets.root"    ,1.41339e+06 , 62140.0),
##    "ZJets_udscg_ptH":("cHists/Hbb_Zmumu__DYJetsToLL_pT100_udscgJets.root"   ,2.64e+06    , 40510.0),
##    "ZJets_bJets_ptH":("cHists/Hbb_Zmumu__DYJetsToLL_pT100_bJets.root"       ,2.64e+06    , 40510.0),
##    "ZZ"             :("cHists/Hbb_Zmumu__ZZ.root"                           ,4.429895e+06, 8255.61),
##    "WW"             :("cHists/Hbb_Zmumu__WW.root"                           ,4.341143e+06, 57109.7),
##    "WZ"             :("cHists/Hbb_Zmumu__WZ.root"                           ,3.764351e+06, 32316.1),
##    "TTbar"          :("cHists/Hbb_Zmumu__TTJets.root"                       ,6.27611e+06,  225197.),
##    "ST-s"           :("cHists/Hbb_Zmumu__T_s-channel.root"                  ,2.39999e+05,  3893.94),
##    "STbar-s"        :("cHists/Hbb_Zmumu__Tbar_s-channel.root"               ,1.39973e+05,  1757.76),
##    ## "ST-t"           :("cHists/Hbb_Zmumu__T_t-channel.root"                  ,xx,           30004.2),
##    "STbar-t"        :("cHists/Hbb_Zmumu__Tbar_t-channel.root"               ,1.455067e+06,  55531.0),
##    "ST-tW"          :("cHists/Hbb_Zmumu__T_tW-channel-DR.root"              ,4.37656e+05,  11177.3),
##    "STbar-tW"       :("cHists/Hbb_Zmumu__Tbar_tW-channel-DR.root"           ,4.63458e+05,  11177.3),
##    }
