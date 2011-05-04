#ifndef __CINT__
#include <iostream>
#include <sstream>
#include <cmath>
#include "Riostream.h"
#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TChain.h"
#endif

vector<TString> inputFiles_;

int readFilelist(std::string filelist){
  int istat=0; //assume success

  std::string CommentLine="#"; // treat lines that begin with "#" as comment lines
  ifstream in;

  string input_line,word,filename;

  cout << "%readFilelist--Reading: " << filelist << endl;
  in.open(filelist.c_str());
  if (!in){
    cerr << "%Could not open file: " << filelist << endl;
    cerr << "%Terminating program" << endl;
    return -1;
  }

  //  while (1) {
  while (getline(in,input_line)){

    if (!in.good()) break;

    int p1 = input_line.find (CommentLine,0); // returns -1 if not found

    if ( p1 < 0){
      istringstream stream(input_line);
      std::vector<string> elements;

      while (stream >> word) {
	elements.push_back(word);
      }

      if (elements.size() >= 1){
	inputFiles_.push_back(elements[0]);
      }else if (elements.size() == 0){
	cout << "%readFilelist--Blank line encountered -- Treating as EOF" << endl;
	return 0;
      }else {
	cout << "%readFilelist--Error parsing input filelist" << endl;
	return -1;
      }
    }
  }
  return istat;

}
