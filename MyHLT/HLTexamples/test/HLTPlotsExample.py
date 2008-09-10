import FWCore.ParameterSet.Config as cms

nevts=200
histofile="histo.root"
inputfile="/store/relval/CMSSW_2_1_2/RelValQCD_Pt_120_170/GEN-SIM-DIGI-RAW-HLTDEBUG-RECO/STARTUP_V5_HF_10TeV_v1/0000/0A6D28A4-6F77-DD11-B9F8-0030487A3232.root"

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
                               MyTrigger = cms.string("HLT_Jet110"),
                               HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets")
                               )
# paths to run
process.p1 = cms.Path(process.plots)

process.schedule = cms.Schedule(process.p1)
