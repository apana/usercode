import FWCore.ParameterSet.Config as cms

############## user-defined parameters ##########################
nevts=100
histofile="l1bits.root"
outfile="l1results.dat"

inputfile = "/store/mc/Summer09/MinBias/GEN-SIM-RECO/MC_31X_V3-v1/0025/C6F7D0C1-E881-DE11-BA55-00E08133CDA0.root"

################################################################

process= cms.Process('L1BitAnalysis')

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevts)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)

# use TFileService for output histograms
process.TFileService = cms.Service("TFileService",
                              fileName = cms.string(histofile)
                              )

process.load("MyL1.L1Analyzer.l1bits_cfi")
process.l1bits.Outfile= cms.string(outfile)


# paths to run
process.p1 = cms.Path(process.l1bits)

process.schedule = cms.Schedule(process.p1)
