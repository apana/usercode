import FWCore.ParameterSet.Config as cms

process = cms.PSet()

# from Analysis.FWLite.files_r179828_L1MuHPF_HLT_unpre import inputFiles
# from Analysis.FWLite.files_r179828_L1EGHPF_HLT_unpre import inputFiles
# from Analysis.FWLite.files_r179828_ZeroBias_PFCorJets import inputFiles
from Analysis.FWLite.files_r179828_ZeroBias_PFCorJets_v2 import inputFiles

myRefTrigger="HLT_ZeroBias_part0_v1"  # reference trigger

myTrigger="HLT_PFJet40_v1"  # trigger
# myTrigger="HLT_PFL2L3CorJet40_v1"
# myTrigger="HLT_PFL1FastL2L3CorJet40_v1"
# myTrigger="HLT_PFL1FastL2L3CorJet40_chs_v1"

Prescl=445  # prescale factor of reference trigger

if (myTrigger == "HLT_PFJet40_v1"):
    WhichJets='hltPFJetsL1Matched'
elif (myTrigger == "HLT_PFL2L3CorJet40_v1"):
    WhichJets='hltPFL2L3CorJetsL1Matched'
elif (myTrigger == "HLT_PFL1FastL2L3CorJet40_v1"):
    WhichJets='hltPFL1FastL2L3CorJetsL1Matched'
elif (myTrigger == "HLT_PFL1FastL2L3CorJet40_chs_v1"):
    WhichJets='hltPFL1FastL2L3CorJetsL1Matchedchs'

# myOutputFile="TriggerRates_" + myTrigger + "_" + WhichJets + "_v2.root"
myOutputFile="Debug.root"

# nevts=100000
nevts=-1

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFiles), ## mandatory
    maxEvents   = cms.int32(nevts),           ## optional
    outputEvery = cms.uint32(10000),       ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('TriggerRates.root'),  ## mandatory
)
process.fwliteOutput.fileName=myOutputFile


process.triggerAnalyzer = cms.PSet(
    ## input specific for this analyzer
    RefTriggerResults = cms.InputTag('TriggerResults::HLT'),
    RefTrigger = cms.string(myRefTrigger),
    Trigger    = cms.string(myTrigger),
    RefLumi    = cms.double(5.37e31),
    TargetLumi = cms.double(5.37e31),
    TriggerResults = cms.InputTag('TriggerResults::TEST'),
    TriggerPS = cms.int32(Prescl), 
    InputJets = cms.InputTag( WhichJets + '::TEST' ),
    BeginLumi=cms.uint32(0),
    EndLumi=  cms.uint32(10000)
)
