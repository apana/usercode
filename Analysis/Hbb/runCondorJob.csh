#!/bin/csh

set noglob   #  prevent filename expansion with wildcard pattern matching

        if ($#argv < 1) then
            echo "Please supply process number"
            echo " "
            exit 1
        endif

# set cfgFile=$argv[1]


set pid=$argv[1]
set workDir=$ANALYZEDIRECTORY
set cfgFile = cfgs/Zee.cfg

set isData=0
set nevts=-1
set doHeavy=0
set doLight=0

set outdir=2012hsts_new

switch ($pid)
    case 0:
        set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DoubleElectron_*.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DoubleElectron2012_all.root
    	set isData = 1
    	breaksw
    case 1:
        set doLight=1
    	set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_M-50_udscgJets.root
    	breaksw
    case 2:
        set doHeavy=1
    	set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Skim_bjets.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_M-50_bJets.root
    	breaksw
    case 3:
        # set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_TTJets_TuneZ2star_8TeV-madgraph-tauola.root"
	set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_TTJets_Merged.root"
    	set outputFiles = "$outdir"/Hbb_Zee__TTJets.root
    	breaksw
    case 4:
        set inputFiles  = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp.root"
    	set outputFiles = "$outdir"/Hbb_Zee__HToBB_M-120.root
    	breaksw
    case 5:
        set doLight=1
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_pT100_udscgJets.root
       	breaksw
    case 6:
        set doHeavy=1
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph_Skim_bjets.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_pT100_bJets.root
       	breaksw
    case 7:
        set doLight=1
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_pT70-100_udscgJets.root
       	breaksw
    case 8:
        set doHeavy=1
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Skim_bjets.root"
    	set outputFiles = "$outdir"/Hbb_Zee__DYJetsToLL_pT70-100_bJets.root
       	breaksw
    case 9:
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_WW_TuneZ2star_8TeV_pythia6_tauola.root"
    	set outputFiles = "$outdir"/Hbb_Zee__WW.root
       	breaksw
    case 10:
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_WZ_TuneZ2star_8TeV_pythia6_tauola.root"
    	set outputFiles = "$outdir"/Hbb_Zee__WZ.root
       	breaksw
    case 11:
    	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_ZZ_TuneZ2star_8TeV_pythia6_tauola.root"
    	set outputFiles = "$outdir"/Hbb_Zee__ZZ.root
       	breaksw
    case 12:
	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root"
	set outputFiles = "$outdir"/Hbb_Zee__T_tW-channel-DR.root
       	breaksw    
    case 13:
	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root"
	set outputFiles = "$outdir"/Hbb_Zee__Tbar_tW-channel-DR.root
        breaksw
    case 14:
	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_T_s-channel_TuneZ2star_8TeV-powheg-tauola.root"
	set outputFiles = "$outdir"/Hbb_Zee__T_s-channel.root
       	breaksw
    case 15:
	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola.root"
	set outputFiles = "$outdir"/Hbb_Zee__Tbar_s-channel.root
       	breaksw
    case 16:
	set inputFiles = "/pnfs/cms/WAX/11/store/user/lpchbb/apana/Step1V33_Step2_V2/DiJetPt_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola.root"
	set outputFiles = "$outdir"/Hbb_Zee__Tbar_t-channel.root
       	breaksw
    default:
        echo "not setup for this"
        exit 1
endsw

echo "Beginning runCondorJob.csh"

echo "-------------------------------"
echo "Current Directory: "
pwd
echo "-------------------------------"

source /uscmst1/prod/sw/cms/setup/cshrc prod
setenv SCRAM_ARCH slc5_amd64_gcc462
cd $workDir

eval `scram runtime -csh`

echo "-------------------------------"
echo "Working Directory: "
pwd
echo "-------------------------------"


echo "Submitting job on `date`" 
echo
echo "Configuration file: $cfgFile"
echo "--------------------------------------------"
cat "$cfgFile"
echo -n "--------------------------------------------"

if ( $isData == 0) then
    if ( $doLight == 1) then
	python HbbAnalyzer.py $cfgFile $outputFiles $inputFiles -n $nevts --lj
    else if ( $doHeavy == 1) then
	python HbbAnalyzer.py $cfgFile $outputFiles $inputFiles -n $nevts --hj
    else
	python HbbAnalyzer.py $cfgFile $outputFiles $inputFiles -n $nevts
    endif
else
    python HbbAnalyzer.py $cfgFile $outputFiles $inputFiles -n $nevts --data
endif

echo "Job finished on `date`" 
