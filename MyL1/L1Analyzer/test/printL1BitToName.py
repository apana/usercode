import FWCore.ParameterSet.Config as cms


GLOBALTAG = 'START3X_V25B::All'

process = cms.Process("L1BitToName")


### Input source ###################################################

# inputfile="/store/mc/Spring10/MinBias/GEN-SIM-RAW/START3X_V25B-v1/0104/FECFDECD-9739-DF11-A00E-001A92971AAA.root"
inputfile="file:L1AlgoSkim_MinimumBiasMC.root"
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(inputfile)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

############# import of standard configurations ####################

process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

#process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
#process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
process.GlobalTag.globaltag = GLOBALTAG

###################################################################

# unpack the GT

import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
process.hltGtDigis = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()
process.hltGtDigis.DaqGtInputTag = "rawDataCollector"
process.hltGtDigis.UnpackBxInEvent = 5

###################################################################

process.l1bittoname = cms.EDAnalyzer("BitNumbertoName",
                                     L1GtRecordInputTag  = cms.InputTag("hltGtDigis"),
                                     L1GtReadoutRecordInputTag = cms.InputTag("hltGtDigis"),
                                     BitsAndPrescales=cms.string("BitsAndPrescales.txt")
                                     )
#process.l1bittoname.L1GtRecordInputTag = "hltGtDigis"
#process.l1bittoname.L1GtReadoutRecordInputTag = "hltGtDigis"

process.p = cms.Path(process.hltGtDigis + process.l1bittoname)


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.debugModules = ['l1bittoname']
process.MessageLogger.cout = cms.untracked.PSet(
    INFO = cms.untracked.PSet(
        limit = cms.untracked.int32(-1)
    ),
    threshold = cms.untracked.string('DEBUG'), ## DEBUG 

    DEBUG = cms.untracked.PSet( ## DEBUG, all messages  

        limit = cms.untracked.int32(-1)
    )
)
