import FWCore.ParameterSet.Config as cms

process = cms.PSet()


# from Analysis.FWLite.files_r180250_Jet import inputFiles
# from Analysis.FWLite.files_r180250_Jet_PFJet import inputFiles
# from Analysis.FWLite.files_r180250_Comm import inputFiles
# from Analysis.FWLite.files_r179828_L1JetHPF_PFCorJets_v2 import inputFiles
from Analysis.FWLite.files_r179828_L1JetHPF_PFCorJets_v3 import inputFiles

# from Analysis.FWLite.files_r179828_ZeroBias_PFCorJets_v2 import inputFiles

# myTrigger="HLT_Jet300_v9"  # reference trigger
# Prescl=1  # prescale factor of reference trigger

# myTrigger="HLT_Jet190_v9"  # reference trigger
# Prescl=100  # prescale factor of reference trigger

# myTrigger="HLT_Jet110_v9"  # reference trigger
# Prescl=1600  # prescale factor of reference trigger

# myTrigger="HLT_Jet60_v9"  # reference trigger
# Prescl=32000  # prescale factor of reference trigger

# myTrigger="HLT_Jet30_v9"  # reference trigger
# Prescl=1440000  # prescale factor of reference trigger

myTrigger="Any"  # reference trigger
# myTrigger="HLT_L1SingleJet92_v4"  # reference trigger
Prescl=1  # prescale factor of reference trigger

# myTrigger="HLT_ZeroBias_part0_v1"  # reference trigger
# Prescl=1  # prescale factor of reference trigger

myOutputFile="Histograms_CompJets_HLTRECO_" + myTrigger + "_L1L2L3Corr.root"
# myOutputFile="Debug.root"

nevts=-1
# nevts=10000
process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFiles), ## mandatory
    maxEvents   = cms.int32(nevts),           ## optional
    outputEvery = cms.uint32(10000),       ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('CompJets.root'),  ## mandatory
)
process.fwliteOutput.fileName=myOutputFile


process.compJets = cms.PSet(
    ## input specific for this analyzer
    ## InputJets1 = cms.InputTag('hltAntiKT5PFJets'),
    # InputJets1 = cms.InputTag('hltPFL2L3CorJetsL1Matched'),
    # InputJets2 = cms.InputTag('hltAntiKT5PFJets'),

    InputJets1 = cms.InputTag('hltPFL1FastL2L3CorJetsL1Matched'),
    InputJets2 = cms.InputTag('ak5PFJetsL1L2L3Residual'),
    TriggerResults = cms.InputTag('TriggerResults::HLT'),
    ## TriggerResults = cms.InputTag('TriggerResults::TEST'),
    Trigger   = cms.string(myTrigger),
    TriggerPS = cms.int32(Prescl), 
    Debug     = cms.bool(False),
    BeginLumi = cms.uint32(7),
    EndLumi   = cms.uint32(1000)
)
