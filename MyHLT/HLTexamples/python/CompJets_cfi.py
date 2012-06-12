import FWCore.ParameterSet.Config as cms

compJets = cms.EDAnalyzer("CompJets",
                       MyTrigger = cms.string("NoTrigger"), ## use the trigger "NoTrigger" for no trigger req.
                       HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
                       PFJetCollection1 = cms.InputTag("ak5PFJetsL1FastL2L3"),
                       PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected"),
                       GenJetCollection = cms.InputTag("ak5GenJets"),
                       Debug=cms.bool(False)
                       )
