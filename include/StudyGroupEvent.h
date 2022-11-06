#pragma once

#include <TString.h>
#include <TH1F.h>
#include <map>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "LEAF/Analyzer/include/RecoEvent.h"

using namespace std;

// Container class for all quantities
class StudyGroupEvent : public RecoEvent{

public:
  // Constructors, destructor
  StudyGroupEvent();
  ~StudyGroupEvent();

  void clear();
  void reset();

};
