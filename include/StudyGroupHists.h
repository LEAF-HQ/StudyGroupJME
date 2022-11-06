#pragma once

#include <TString.h>
#include <TH1F.h>
#include <map>
#include <TTreeReader.h>
#include <TROOT.h>
#include "LEAF/StudyGroup/include/StudyGroupEvent.h"
#include "LEAF/Analyzer/include/BaseHists.h"

using namespace std;

class StudyGroupHists : public BaseHists{

public:
  // Constructors, destructor
  StudyGroupHists(TString dir_);
  StudyGroupHists(const StudyGroupHists &) = default;
  StudyGroupHists & operator = (const StudyGroupHists &) = default;
  ~StudyGroupHists() = default;

  // Main functions
  void fill(const StudyGroupEvent & event);


protected:

  shared_ptr<TH1D> hmetpt, hmetphi, hsumweights;

};
