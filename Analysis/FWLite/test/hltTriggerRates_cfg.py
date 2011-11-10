import FWCore.ParameterSet.Config as cms

process = cms.PSet()

# from Analysis.FWLite.files_r179828_L1MuHPF_HLT_unpre import inputFiles
# from Analysis.FWLite.files_r179828_L1EGHPF_HLT_unpre import inputFiles
from Analysis.FWLite.files_r179828_L1JetHPF_HLT_unpre import inputFiles

myRefTrigger="HLT_L1HTT100_v1"  # reference trigger
Prescl=20  # prescale factor of reference trigger

myOutputFile="TriggerRates_" + myRefTrigger + ".root"

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(inputFiles), ## mandatory
    maxEvents   = cms.int32(-1),           ## optional
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
    RefLumi    = cms.double(5.37e31),
    TargetLumi = cms.double(5.37e31),
    TriggerResults = cms.InputTag('TriggerResults::TEST'),
    TriggerPS = cms.int32(Prescl), 
    BeginLumi=cms.uint32(0),
    EndLumi=  cms.uint32(10000)
)
