#include "LEAF/StudyGroup/include/StudyGroupEvent.h"
#include "LEAF/Analyzer/include/constants.h"
#include <TH1D.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TKey.h>
#include <TTree.h>
#include <TLatex.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <iostream>
#include <sys/stat.h>

using namespace std;

StudyGroupEvent::StudyGroupEvent(){
}

StudyGroupEvent::~StudyGroupEvent(){
}

void StudyGroupEvent::clear(){
  RecoEvent::clear();
}

void StudyGroupEvent::reset(){
  RecoEvent::reset();
}
