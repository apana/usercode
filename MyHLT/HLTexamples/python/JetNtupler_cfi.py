import FWCore.ParameterSet.Config as cms

jetNtupler = cms.EDAnalyzer("JetNtupler",
                       HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
                       RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
                       HLTJets  = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected"),
                       GenJets  = cms.InputTag("ak5GenJets"),
                       Rho      = cms.InputTag("kt6PFJets:rho:RECO" ),
                       Vertices = cms.InputTag('offlinePrimaryVertices'),
                       Debug=cms.bool(False),
                       Monte=cms.bool(False)
                       )
