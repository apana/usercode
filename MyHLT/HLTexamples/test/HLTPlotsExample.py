import FWCore.ParameterSet.Config as cms

nevts=-1
histofile="histo.root"
inputfile="/store/relval/2008/4/28/RelVal-RelValQCD_Pt_80_120-1209247429-IDEAL_V1-2nd/0001/0C5D7E6C-0315-DD11-B1F0-000423D99BF2.root"

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
                               MyTrigger = cms.string("HLT1jet80"),
                               HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets")
                               )
# paths to run
process.p1 = cms.Path(process.plots)

process.schedule = cms.Schedule(process.p1)
