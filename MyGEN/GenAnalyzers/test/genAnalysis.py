import FWCore.ParameterSet.Config as cms

process= cms.Process('ANALYSIS')

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/apana/HLTout_EWK_Ztautau_GEN.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.genAnalysis = cms.EDAnalyzer("GenEventAnalyzer",
                                     mctruth = cms.InputTag("genParticles"),
                                     genEventScale = cms.InputTag("genEventScale"),
                                     Histogram = cms.string("genanalysis.root")
                                     )

process.p1 = cms.Path(process.genAnalysis)
process.schedule = cms.Schedule(process.p1)
