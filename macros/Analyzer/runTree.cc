void runTree(){


  gROOT->ProcessLine(".L processTree.cc++");

  int nevts=-1; // -1 == process all events on tree
  //nevts=5;

  ostringstream ch_n("");
  ch_n << nevts;


  bool doXSweighting=false; 
  bool skimmedFiles=false;
  //skimmedFiles=true;

  string filelist="files_RelVal_160_wisc.list";
  string histname="hlt_hists_RelVal_160";
  if (nevts>0) histname = histname + "_n" + ch_n.str();
  if (! doXSweighting) histname = histname + "_unwt";
  histname = histname + ".root";

  // run the job
  processTree(filelist,histname,nevts,doXSweighting,skimmedFiles);

}
