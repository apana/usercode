import FWCore.ParameterSet.Config as cms

nevts=200
histofile="l1bits.root"
# inputfile="/store/data/Commissioning08/Calo/RECO/v1/000/068/288/92E74EA6-63A8-DD11-8628-001617DBCF1E.root"
# inputfile="file:/uscmst1b_scratch/lpc1/lpctrig/apana/dev/CMSSW_3_1_0_pre9/src/run/RelVal_HLT2_8E29_HF07_v1_HF07_500ev_QCD.root"
inputfile="/store/relval/CMSSW_3_1_0/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/STARTUP31X_V1-v1/0001/F472473A-9966-DE11-A190-000423D6A6F4.root"

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
                               Outfile   = cms.string("l1results_8e29.dat"),
                               # L1GtRecordInputTag  = cms.InputTag("gtDigis"),
                               L1GtRecordInputTag  = cms.InputTag("hltGtDigis::HLT"),
                               L1GtObjectMapRecord = cms.InputTag("hltL1GtObjectMap::HLT"),
                               L1AlgoName = cms.string("L1_SingleJet6"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets")
                               )
# paths to run
process.p1 = cms.Path(process.l1plots)

process.schedule = cms.Schedule(process.p1)
