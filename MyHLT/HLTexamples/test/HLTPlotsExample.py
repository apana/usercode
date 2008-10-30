import FWCore.ParameterSet.Config as cms

nevts=200
histofile="histo.root"
inputfile="file:/uscmst1b_scratch/lpc1/lpctrig/apana/fsim/CMSSW_2_1_10/src/run/AODIntegrationTestWithHLT.root"

process= cms.Process('HLTPlots')

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevts)
)

# use TFileService for output histograms
process.TFileService = cms.Service("TFileService",
                              fileName = cms.string(histofile)
                              )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)

process.plots = cms.EDAnalyzer("HLTPlotsExample",
                               # MyTrigger = cms.string("HLT_Jet110"),
                               MyTrigger = cms.string("HLT_L1MuOpen"),
                               # MyTrigger = cms.string("HLT_Mu3"),
                               HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets"),
                               MuonCollection   = cms.InputTag("muons")
                               )
# paths to run
process.p1 = cms.Path(process.plots)

process.schedule = cms.Schedule(process.p1)
