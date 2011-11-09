import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_10_1_MiA.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_11_1_Jn0.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_12_1_ANx.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_13_1_1mM.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_14_1_9Ce.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_15_1_f92.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_16_1_RA2.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_17_1_TPb.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_18_1_2fD.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_19_1_cfZ.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_1_1_cbg.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_20_1_mWk.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_21_1_fZ4.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_22_1_Mzq.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_23_1_FuQ.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_24_1_YmX.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_25_1_kiu.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_26_1_oSQ.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_27_1_5ks.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_28_1_tdq.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_2_1_QS7.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_3_1_YVJ.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_4_1_WPA.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_5_1_kdr.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_6_1_eUR.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_7_1_dvv.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_8_1_rc9.root",
        "dcap://cmsdca.fnal.gov:24136/pnfs/fnal.gov/usr/cms/WAX/11/store/user/apana/r179828_L1JetHPF/L1AlgoSkim_3e33_9_1_Q58.root"
        ), ## mandatory
    maxEvents   = cms.int32(100),                             ## optional
    outputEvery = cms.uint32(1000),                            ## optional
)
    
process.fwliteOutput = cms.PSet(
    fileName  = cms.string('TriggerRates.root'),  ## mandatory
)

process.triggerAnalyzer = cms.PSet(
    ## input specific for this analyzer
    RefTriggerResults = cms.InputTag('TriggerResults::HLT'),
    RefTrigger = cms.string('HLT_L1SingleJet92_v4'),
    TriggerResults = cms.InputTag('TriggerResults::TEST'),
    TriggerPS = cms.int32(1),
    BeginLumi=cms.uint32(0),
    EndLumi=  cms.uint32(10000)
)
