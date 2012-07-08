import sys,string
import FWCore.ParameterSet.Config as cms

######################################################################

nevts=-1
# nevts=100

DS="Jet"
## DS="Commish"
## DS="RelValQCD"
## DS="QCD"

WhichPreTrigger='HLT_*'  # accept all HLT paths
## WhichPreTrigger='HLT_L1SingleJet16_v*'
## WhichPreTrigger='HLT_PFJet260_v*'

# WhichJEC="GT"
# WhichJEC="Summer12_V7_MC"
# WhichJEC="Summer12_V7_DATA"
WhichJEC="Jec12_V7"

histofile="JetNtupler_" + DS
if WhichPreTrigger != "HLT_*":
    histofile=histofile + "_Pre" + WhichPreTrigger[:WhichPreTrigger.rfind("_v")]
histofile=histofile + "_r195655_HLTJECv8_OfflineJECfrom"+WhichJEC+"_new.root"
# histofile="tmp.root"

if DS == "Jet":
    # from MyHLT.HLTexamples.files_r191226_Jet import inputFiles
    from MyHLT.HLTexamples.files_r195655_JetDS_v8 import inputFiles
elif DS == "Commish":
    # from MyHLT.HLTexamples.files_r191226_Commish import inputFiles
    from MyHLT.HLTexamples.files_r195655_CommishDS_v8 import inputFiles
elif DS == "RelValQCD":
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000 import inputFiles
    from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_PU_JECv7 import inputFiles
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_noPU import inputFiles
elif DS == "QCD":
    from MyHLT.HLTexamples.files_QCD_FlatPt_15to3000_lungu import inputFiles
else:
    print "Bad DS: ",DS
    sys.exit(1)

# inputFiles=["file:../../../run/outputHLTDQMResults_0.3.root"]
# histofile="JetNtupler_03.root"

################################################################################

process= cms.Process('JetNtupler')

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nevts)
)

# use TFileService for output histograms
process.TFileService = cms.Service("TFileService",
                              fileName = cms.string(histofile)
                              )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.GlobalTag.globaltag = "GR_R_52_V9::All"

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles)
)

#####################################################
##### setup several instances of the JetNtupler analyzer
######################################################

from MyHLT.HLTexamples.JetNtupler_cfi import jetNtupler
# ccla
process.jetNtupler_uncorr = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJets"),
    # Rho      = cms.InputTag("hltKT6PFJets:rho:TEST" ),
    HLTJets = cms.InputTag("hltAntiKT5PFJets")
)

process.jetNtupler_uncorrNoPUHLT_L1L2L3RECO = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
    HLTJets = cms.InputTag("hltAntiKT5PFJetsNoPU")
)

process.jetNtupler_L1L2L3chsHLT_L1L2L3RECO = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU")
)

process.jetNtupler_uncorrHLT_L1L2L3RECO = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
    HLTJets = cms.InputTag("hltAntiKT5PFJets")
)

process.jetNtupler_L1L2L3 = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected")
)

process.jetNtupler_L1L2L3chs = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFchsJetsL1FastL2L3"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU")
)


process.jetNtupler_L1L2L3ResidchsReco_L1L2L3chsHLT = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFchsJetsL1FastL2L3Residual"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU")
)

process.jetNtupler_L1L2L3ResidchsReco_L1L2L3HLT = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFchsJetsL1FastL2L3Residual"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected")
)

process.jetNtupler_L1L2L3ResidReco_L1L2L3HLT = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3Residual"),
    HLTJets = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected")
)


######################################################################
### Run jet corrections for RECO jets
######################################################################

process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True


###
if WhichJEC.find("Summer")==0:

    SQLFile="sqlite_file:"+WhichJEC +".db"
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                           connect = cms.string(SQLFile),
                           # timetype = cms.string('runnumber'),
                           # cms.string("frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS"),
                           toGet =  cms.VPSet(
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_" + WhichJEC + "_AK5Calo"),
                     label= cms.untracked.string("AK5Calo")),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_" + WhichJEC + "_AK5PF"),
                     label=cms.untracked.string("AK5PF")),
            cms.PSet(record = cms.string("JetCorrectionsRecord"),
                     tag = cms.string("JetCorrectorParametersCollection_" + WhichJEC + "_AK5PFchs"),
                     label=cms.untracked.string("AK5PFchs")),
            
            )
                               )
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


process.pfPileUp = cms.EDProducer("PFPileUp",
    PFCandidates = cms.InputTag("particleFlow"),
    Enable = cms.bool(True),
    checkClosestZVertex = cms.bool(False),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("offlinePrimaryVertices")
)

process.pfNoPileUp = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlow"),
    enable = cms.bool(True),
    topCollection = cms.InputTag("pfPileUp"),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    verbose = cms.untracked.bool(False)
)

process.ak5PFchsJets = cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(True),
    voronoiRfact = cms.double(-0.9),
    maxBadHcalCells = cms.uint32(9999999),
    doAreaDiskApprox = cms.bool(False),
    maxRecoveredEcalCells = cms.uint32(9999999),
    jetType = cms.string('PFJet'),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    doRhoFastjet = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    maxBadEcalCells = cms.uint32(9999999),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    maxRecoveredHcalCells = cms.uint32(9999999),
    rParam = cms.double(0.5),
    maxProblematicHcalCells = cms.uint32(9999999),
    doOutputJets = cms.bool(True),
    src = cms.InputTag("pfNoPileUp"),
    inputEtMin = cms.double(0.0),
    srcPVs = cms.InputTag(""),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    maxProblematicEcalCells = cms.uint32(9999999),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0)
)

process.makeAK5PFNoPUJets = cms.Sequence(process.pfPileUp + process.pfNoPileUp + process.ak5PFchsJets)

process.ak5PFchsL1Fastjet = cms.ESProducer( "L1FastjetCorrectionESProducer",
    era         = cms.string('Summer11'),
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK5PFchs'),
    section     = cms.string(''),
    srcRho      = cms.InputTag('kt6PFJets','rho'),
    useCondDB = cms.untracked.bool(True)
    )

process.ak5PFchsL2L3 = cms.ESProducer( "JetCorrectionESChain",
    correctors = cms.vstring('ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )
process.ak5PFchsL1FastL2L3 = cms.ESProducer( "JetCorrectionESChain",
    correctors = cms.vstring('ak5PFchsL1Fastjet','ak5PFchsL2Relative','ak5PFchsL3Absolute')
    )


process.ak5PFchsL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFchsL1Fastjet', 
        'ak5PFchsL2Relative', 
        'ak5PFchsL3Absolute', 
        'ak5PFchsResidual')
)

process.ak5PFchsL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2Relative')
)


process.ak5PFchsL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L3Absolute')
)

process.ak5PFchsResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2L3Residual')
)

process.ak5PFchsJetsL1FastL2L3 = cms.EDProducer(
    'PFJetCorrectionProducer',
    src         = cms.InputTag('ak5PFchsJets'),
    correctors  = cms.vstring('ak5PFchsL1FastL2L3')
    ) 

process.ak5PFchsJetsL1FastL2L3Residual = cms.EDProducer(
    'PFJetCorrectionProducer',
    src         = cms.InputTag('ak5PFchsJets'),
    correctors  = cms.vstring('ak5PFchsL1FastL2L3Residual')
    ) 


###
######################################################################
### import triggerResultsFilter to preFilter on HLT paths if desired
######################################################################

from HLTrigger.HLTfilters.triggerResultsFilter_cfi import triggerResultsFilter

process.trigPreFilter = triggerResultsFilter.clone(
    hltResults              = cms.InputTag('TriggerResults::HLT'),
    triggerConditions       = cms.vstring(WhichPreTrigger)
    )


process.p1 = cms.Path(process.trigPreFilter *
                      process.ak5PFJetsL2L3 *
                      process.ak5PFJetsL1FastL2L3 *
                      process.ak5PFJetsL1FastL2L3Residual * 
                      process.makeAK5PFNoPUJets *
                      process.ak5PFchsJetsL1FastL2L3 *
                      process.ak5PFchsJetsL1FastL2L3Residual *
                      process.jetNtupler_uncorr *
                      process.jetNtupler_uncorrHLT_L1L2L3RECO *
                      process.jetNtupler_uncorrNoPUHLT_L1L2L3RECO *
                      process.jetNtupler_L1L2L3chsHLT_L1L2L3RECO *
                      process.jetNtupler_L1L2L3 *
                      process.jetNtupler_L1L2L3chs *
                      process.jetNtupler_L1L2L3ResidReco_L1L2L3HLT *
                      process.jetNtupler_L1L2L3ResidchsReco_L1L2L3HLT *
                      process.jetNtupler_L1L2L3ResidchsReco_L1L2L3chsHLT
                      )

print "\n\tHistograms written to: ",histofile,"\n"

## temp = process.dumpPython()
## outputfile = file("expstd.py",'w')
## outputfile.write(temp)
## outputfile.close()
