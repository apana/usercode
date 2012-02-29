import FWCore.ParameterSet.Config as cms

process = cms.PSet()


# from Analysis.FWLite.files_r180250_Jet import inputFiles
from Analysis.FWLite.files_r180250_Jet_PFJet import inputFiles
# from Analysis.FWLite.files_r180250_Comm import inputFiles

# myTrigger="HLT_Jet300_v9"  # reference trigger
# Prescl=1  # prescale factor of reference trigger

# myTrigger="HLT_Jet190_v9"  # reference trigger
# Prescl=100  # prescale factor of reference trigger

# myTrigger="HLT_Jet110_v9"  # reference trigger
# Prescl=1600  # prescale factor of reference trigger

# myTrigger="HLT_Jet60_v9"  # reference trigger
# Prescl=32000  # prescale factor of reference trigger

myTrigger="HLT_Jet30_v9"  # reference trigger
Prescl=1440000  # prescale factor of reference trigger

# myTrigger="HLT_L1SingleJet16_v4"  # reference trigger
# Prescl=3600000  # prescale factor of reference trigger


myOutputFile="Histograms_CompJets_" + myTrigger + "_L2L2Corr.root"

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFiles), ## mandatory
    maxEvents   = cms.int32(-1),           ## optional
    outputEvery = cms.uint32(10000),       ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('CompJets.root'),  ## mandatory
)
process.fwliteOutput.fileName=myOutputFile


process.compJets = cms.PSet(
    ## input specific for this analyzer
    ## InputJets1 = cms.InputTag('hltAntiKT5PFJets'),
    InputJets1 = cms.InputTag('ak5PFJetsL2L3Residual'),
    InputJets2 = cms.InputTag('xxx'),
    TriggerResults = cms.InputTag('TriggerResults::HLT'),
    ## TriggerResults = cms.InputTag('TriggerResults::TEST'),
    Trigger = cms.string(myTrigger),
    TriggerPS = cms.int32(Prescl), 
    BeginLumi=cms.uint32(7),
    EndLumi=  cms.uint32(1000)
)
