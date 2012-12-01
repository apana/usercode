import FWCore.ParameterSet.Config as cms

ptcut=300

InputFile="dcache:/pnfs/cms/WAX/11/store/user/lpchbb/apana/dev/Skim_pT300/ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp_Skim_pTH300_3_1_O4n.root"

OutputFile="Skim_Higgs"+str(ptcut)+".root"

nevts=20



process = cms.Process("HLTSkim")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(InputFile),
    skipBadFiles = cms.untracked.bool(True)                  
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( nevts )
)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
import MyGEN.SkimGen.SkimGen_cfi as MCSkim
#

process.mySkim = MCSkim.skim.clone(
    # genparts=cms.vint32(25), # Higgs
    genparts=cms.vint32(23,24), # W,Z
    genpt=cms.double(ptcut),
    Debug=cms.bool(False)
    )

process.FilterPath = cms.Path( process.mySkim )

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'FilterPath') ), 
    fileName = cms.untracked.string(OutputFile)
)

process.o = cms.EndPath( process.out )
# process.schedule = cms.Schedule(process.o)
