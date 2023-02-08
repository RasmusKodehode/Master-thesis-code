#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <TTree.h>
#include <vector>
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <TH1.h>
#include <AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h>
#include <AsgTools/ToolHandle.h>

class MyxAODAnalysis : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);
    ToolHandle<IGoodRunsListSelectionTool> m_grl;

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  /// output variables for the current event
  /// \{
  unsigned int m_runNumber = 0; ///< Run number
  unsigned long long m_eventNumber = 0; ///< Event number
  /// Muon 4-momentum variables
  std::vector<float> *m_muonEta = nullptr;
  std::vector<float> *m_muonPhi = nullptr;
  std::vector<float> *m_muonPt = nullptr;
  std::vector<float> *m_muonE = nullptr;
  std::vector<float> *m_leadingPt = nullptr;
  std::vector<float> *m_subPt = nullptr;
  std::vector<float> *m_muonSize = nullptr;
  std::vector<float> *m_Z_mass = nullptr;
  std::vector<float> *m_truthZ = nullptr;
  std::vector<float> *m_children = nullptr;
  std::vector<float> *m_mChildren = nullptr;
  std::vector<float> *m_Z_match = nullptr;
  std::vector<float> *m_truthEta = nullptr;
  std::vector<float> *m_truthPhi = nullptr;
  std::vector<float> *m_matchEta = nullptr;
  std::vector<float> *m_matchPhi = nullptr;
  std::vector<float> *m_ptRes = nullptr;
  std::vector<float> *m_etaRes = nullptr;
  std::vector<float> *m_MinBias = nullptr;
  /// \}

private:
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;
};

#endif