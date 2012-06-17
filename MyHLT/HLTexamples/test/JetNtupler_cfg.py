import sys,string
import FWCore.ParameterSet.Config as cms

######################################################################

nevts=-10000
# nevts=100

DS="Jet"
## DS="Commish"
## DS="RelValQCD"
## DS="QCD"

## WhichPreTrigger='HLT_*'  # accept all HLT paths
## WhichPreTrigger='HLT_L1SingleJet16_v*'
WhichPreTrigger='HLT_PFJet260_v*'


histofile="JetNtupler_" + DS
if WhichPreTrigger != "HLT_*":
    histofile=histofile + "_Pre" + WhichPreTrigger[:WhichPreTrigger.rfind("_v")]
histofile=histofile + "_r1933366_JECv7.root"

if DS == "Jet":
    # from MyHLT.HLTexamples.files_r191226_Jet import inputFiles
    from MyHLT.HLTexamples.files_r193336_Jet import inputFiles
elif DS == "Commish":
    # from MyHLT.HLTexamples.files_r191226_Commish import inputFiles
    from MyHLT.HLTexamples.files_r193336_Commish import inputFiles
elif DS == "RelValQCD":
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000 import inputFiles
    from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_PU_JECv7 import inputFiles
    # from MyHLT.HLTexamples.files_RelValQCD_FlatPt_15_3000_noPU import inputFiles
elif DS == "QCD":
    from MyHLT.HLTexamples.files_QCD_FlatPt_15to3000_lungu import inputFiles
else:
    print "Bad DS: ",DS
    sys.exit(1)


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
process.GlobalTag.globaltag = "GR_R_52_V7::All"

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputFiles)
)

#####################################################
##### setup several instances of the JetNtupler analyzer
######################################################

from MyHLT.HLTexamples.JetNtupler_cfi import jetNtupler

process.jetNtupler_uncorr = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJets"),
    HLTJets = cms.InputTag("hltAntiKT5PFJets")
)

process.jetNtupler_corr = jetNtupler.clone(
    HLTriggerResults = cms.InputTag("TriggerResults::HLT"),
    RecoJets = cms.InputTag("ak5PFJetsL1FastL2L3"),
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
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                           connect = cms.string('sqlite_file:Jec12_V7.db'),
                           # timetype = cms.string('runnumber'),
                           # cms.string("frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS"),
                           toGet =  cms.VPSet(
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5Calo"),
             label= cms.untracked.string("AK5Calo")),
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5PF"),
             label=cms.untracked.string("AK5PF")),
    cms.PSet(record = cms.string("JetCorrectionsRecord"),
             tag = cms.string("JetCorrectorParametersCollection_Jec12_V7_AK5PFchs"),
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
                      # process.ak5PFJetsL1FastL2L3Residual * 
                      process.jetNtupler_uncorr *
                      process.jetNtupler_corr
                      )

print "\n\tHistograms written to: ",histofile,"\n"

## temp = process.dumpPython()
## outputfile = file("expstd.py",'w')
## outputfile.write(temp)
## outputfile.close()
