import sys,string
import FWCore.ParameterSet.Config as cms

######################################################################


OutputFile="Jets_SQL.root"
nevts=-1000
# nevts=10000

DS="Jet"
## DS="Commish"
## DS="RelValQCD"
## DS="QCD"

WhichPreTrigger='HLT_*'  # accept all HLT paths
## WhichPreTrigger='HLT_L1SingleJet16_v*'

# set to NoTrigger for no trigger req in CompJets
WhichTrigger="HLT_PFJet320_v5"  
## WhichTrigger="NoTrigger"


histofile="compJets_" + WhichTrigger + "_" + DS
if WhichPreTrigger != "HLT_*":
    histofile=histofile + "_Pre" + WhichPreTrigger[:WhichPreTrigger.rfind("_v")]
# histofile=histofile + "_PU_newJEC.root"
# histofile=histofile + "_PU_JECv7v7.root"
# histofile=histofile + "_r191226_JECv7_tst.root"
histofile=histofile + "_r193336_JECv8_OfflinefromMCSQL.root"

if DS == "Jet":
    # from MyHLT.HLTexamples.files_r191226_Jet import inputFiles
    from MyHLT.HLTexamples.files_r193336_JetDS_v8 import inputFiles
elif DS == "Commish":
    # from MyHLT.HLTexamples.files_r191226_Commish import inputFiles
    from MyHLT.HLTexamples.files_r193336_CommishDS_v8 import inputFiles
elif DS == "RelValQCD":
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000 import inputFiles
    from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_PU_JECv7 import inputFiles
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_noPU import inputFiles
elif DS == "QCD":
    from MyHLT.HLTexamples.files_QCD_FlatPt_15to3000_lungu import inputFiles
else:
    print "Bad DS: ",DS
    sys.exit(1)

# inputFiles=["file:../../../run/outputHLTDQMResults.root"]
# inputFiles=['dcache:/pnfs/cms/WAX/11/store/user/apana/PFCorJets/r193336_chs_JEC_v8_Jet/outputHLTDQMResults_5_1_0TT.root']

# histofile="tmp.root"
################################################################################

process= cms.Process('CompJets')

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
process.GlobalTag.globaltag = "GR_R_52_V7::All"
# process.GlobalTag.globaltag = "GR_R_52_V9::All"

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles)
)

#####################################################
##### setup several instances of the CompJets analyzer
######################################################

from MyHLT.HLTexamples.CompJets_cfi import compJets

process.compJets_uncorr = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFJets"),
    PFJetCollection2 = cms.InputTag("hltAntiKT5PFJets")
)

process.compJets_corr = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFJetsL1FastL2L3"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected"),
    Debug            = cms.bool(False)
)

process.compJets_chs_corr = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFchsJetsL1FastL2L3"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU"),
    Debug            = cms.bool(False)
)

process.compJets_chs_corrWresid = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFchsJetsL1FastL2L3Residual"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU"),
    Debug            = cms.bool(False)
)

process.compJets_corrWresid = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFJetsL1FastL2L3Residual"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected")
)

process.compJets_corrL2L3reco = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFJetsL2L3"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL2L3Corrected")
)

process.compJets_recoCHS = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("ak5PFJetsL1FastL2L3"),
    PFJetCollection2 = cms.InputTag("ak5PFchsJetsL1FastL2L3")
)

process.compJets_hltCHS = compJets.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    PFJetCollection1 = cms.InputTag("hltAK5PFJetL1FastL2L3Corrected"),
    PFJetCollection2 = cms.InputTag("hltAK5PFJetL1FastL2L3CorrectedNoPU")
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

##

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

## process.ak5PFchsJetsL1FastL2L3 = process.ak5PFchsJetsL2L3.clone(src = 'ak5PFchsJets', correctors = ['ak5PFchsL1FastL2L3'])

## 
###
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                           # connect = cms.string('sqlite_file:Jec12_V7.db'),
                           connect = cms.string('sqlite_file:Summer12_V7_MC.db'),
                           # timetype = cms.string('runnumber'),
                           # cms.string("frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS"),
                           toGet =  cms.VPSet(
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             # tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5Calo"),
             tag = cms.string("JetCorrectorParametersCollection_Summer12_V7_MC_AK5Calo"),
             label= cms.untracked.string("AK5Calo")),
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             # tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5PF"),
             tag = cms.string("JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF"),
             label=cms.untracked.string("AK5PF")),
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             # tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5PFchs"),
             tag = cms.string("JetCorrectorParametersCollection_Summer12_V7_MC_AK5PFchs"),
             label=cms.untracked.string("AK5PFchs")),
           
    )
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


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
                      process.compJets_uncorr *
                      process.compJets_corr *
                      process.compJets_corrL2L3reco *
                      process.compJets_corrWresid *
                      process.compJets_chs_corr *
                      process.compJets_chs_corrWresid *
                      process.compJets_recoCHS *
                      process.compJets_hltCHS
                      )

print "\n\tHistograms written to: ",histofile,"\n"

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    # SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'FilterPath') ), 
    fileName = cms.untracked.string(OutputFile)
)

# process.o = cms.EndPath( process.out )

temp = process.dumpPython()
outputfile = file("expstd.py",'w')
outputfile.write(temp)
outputfile.close()
