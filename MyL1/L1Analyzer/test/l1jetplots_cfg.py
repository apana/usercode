import FWCore.ParameterSet.Config as cms

nevts=200
histofile="l1jetplots.root"
# inputfile="/store/data/Commissioning08/Calo/RECO/v1/000/068/288/92E74EA6-63A8-DD11-8628-001617DBCF1E.root"
# inputfile="file:/uscmst1b_scratch/lpc1/lpctrig/apana/dev/CMSSW_3_1_0_pre9/src/run/RelVal_HLT2_8E29_HF07_v1_HF07_500ev_QCD.root"
inputfile="/store/relval/CMSSW_3_1_0/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V1-v1/0001/F472473A-9966-DE11-A190-000423D6A6F4.root"

process= cms.Process('L1JetsAnalysis')

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

process.l1jets = cms.EDAnalyzer("L1JetPlots",
                                 CaloJetAlgorithm = cms.InputTag("hltIterativeCone5CaloJets::HLT"),
                                 GenJetAlgorithm  = cms.InputTag("iterativeCone5GenJets"),
                                 recmet           = cms.InputTag("hltMet::HLT"),
                                 genmet           = cms.InputTag("genMetCalo"),
                                 l1collections    = cms.InputTag("hltL1extraParticles::HLT")
                               )
# paths to run
process.p1 = cms.Path(process.l1jets)

process.schedule = cms.Schedule(process.p1)
