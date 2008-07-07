import FWCore.ParameterSet.Config as cms

nevts=-1
histofile="histogram.root"
inputfile="file:/tmp/apana/HLTout_EWK_Wtaunu_GEN.root"
process= cms.Process('ANALYSIS')

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevts)
)
process.genAnalysis = cms.EDAnalyzer("GenEventAnalyzer",
                                     mctruth = cms.InputTag("genParticles"),
                                     genEventScale = cms.InputTag("genEventScale")
                                     # Histogram = cms.string("genanalysis.root")
                                     )

process.TFileService = cms.Service("TFileService",
                              fileName = cms.string(histofile)
                              )

process.p1 = cms.Path(process.genAnalysis)
process.schedule = cms.Schedule(process.p1)
