
import FWCore.ParameterSet.Config as cms

nevts=-1
histofile="histo.root"
# inputfile="/store/data/BeamCommissioning09/MinimumBias/RAW/v1/000/123/596/2EBD6495-39E2-DE11-A0B3-000423D99394.root"
inputfile="/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0034/FC56D7F4-12D2-E011-A1C1-0026189438D5.root"

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
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0034/FC56D7F4-12D2-E011-A1C1-0026189438D5.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0033/42FC4057-8AD0-E011-BD51-00261894397A.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0000/E88DA4F5-73D0-E011-8866-00261894390E.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0000/A2198EF8-73D0-E011-8B1A-002618943918.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0000/909C7E00-74D0-E011-9C9B-002618943937.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0000/883016E7-73D0-E011-A72E-0026189438B3.root",
        "/store/relval/CMSSW_4_2_9_HLT1_hltpatch1/Jet/RECO/GR_R_42_V14_RelVal_jet2010B-v1/0000/3A5A9DED-73D0-E011-BB36-002618943953.root"
        )
)

process.plots = cms.EDAnalyzer("HLTPlotsExample",
                               MyTrigger = cms.string("HLT_Jet15U_v3"),
                               # MyTrigger = cms.string("HLT_L1MuOpen"),
                               # MyTrigger = cms.string("HLT_Mu3"),
                               HLTriggerResults1 = cms.InputTag("TriggerResults","","HLT"),
                               HLTriggerResults2 = cms.InputTag("TriggerResults::reRECO"),
                               CaloJetAlgorithm = cms.InputTag("iterativeCone5CaloJets"),
                               MuonCollection   = cms.InputTag("muons")
                               )
# paths to run
process.p1 = cms.Path(process.plots)

process.schedule = cms.Schedule(process.p1)
