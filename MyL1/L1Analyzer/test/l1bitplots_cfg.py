import FWCore.ParameterSet.Config as cms

nevts=200
histofile="l1bits.root"
inputfile="/store/data/Commissioning08/Calo/RECO/v1/000/068/288/92E74EA6-63A8-DD11-8628-001617DBCF1E.root"

process= cms.Process('L1BitAnalysis')

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

process.l1plots = cms.EDAnalyzer("L1Bits",
                               Outfile   = cms.string("l1results.dat"),
                               L1GtRecordInputTag  = cms.InputTag("gtDigis"),
                               L1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap"),
                               L1AlgoName = cms.string("L1_SingleJet50_noEE"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets")
                               )
# paths to run
process.p1 = cms.Path(process.l1plots)

process.schedule = cms.Schedule(process.p1)
