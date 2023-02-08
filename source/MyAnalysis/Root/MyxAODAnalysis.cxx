#include <AsgTools/MessageCheck.h>
#include <xAODEventInfo/EventInfo.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODTruth/TruthEventContainer.h>
#include <xAODTruth/TruthParticle.h>
#include <MyAnalysis/MyxAODAnalysis.h>

MyxAODAnalysis :: MyxAODAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
  : EL::AnaAlgorithm (name, pSvcLocator),
    m_grl ("GoodRunsListSelectionTool/grl", this)
{
  // declare the tool handle as a property on the algorithm
  declareProperty ("grlTool", m_grl, "the GRL tool");

  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  This is also where you
  // declare all properties for your algorithm.  Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.
}


MyxAODAnalysis :: ~MyxAODAnalysis () {
  delete m_muonEta;
  delete m_muonPhi;
  delete m_muonPt;
  delete m_muonE;
  delete m_leadingPt;
  delete m_subPt;
  delete m_muonSize;
  delete m_Z_mass;
  delete m_truthZ;
  delete m_children;
  delete m_mChildren;
  delete m_Z_match;
  delete m_truthEta;
  delete m_truthPhi;
  delete m_matchEta;
  delete m_matchPhi;
  delete m_ptRes;
  delete m_etaRes;
  delete m_MinBias;
}

StatusCode MyxAODAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  ANA_MSG_INFO ("in initialize");
  ANA_CHECK (m_grl.retrieve());
  ANA_CHECK (book (TTree ("analysis", "My analysis ntuple")));
  TTree* mytree = tree ("analysis");               //define tree and add branches for each variable
  mytree->Branch ("RunNumber", &m_runNumber);
  mytree->Branch ("EventNumber", &m_eventNumber);
  m_muonEta = new std::vector<float>();
  mytree->Branch ("MuonEta", &m_muonEta);
  m_muonPhi = new std::vector<float>();
  mytree->Branch ("MuonPhi", &m_muonPhi);
  m_muonPt = new std::vector<float>();
  mytree->Branch ("MuonPt", &m_muonPt);
  m_muonE = new std::vector<float>();
  mytree->Branch ("MuonE", &m_muonE);
  m_leadingPt = new std::vector<float>();
  mytree->Branch ("leadingPt", &m_leadingPt);
  m_subPt = new std::vector<float>();
  mytree->Branch ("subPt", &m_subPt);
  m_muonSize = new std::vector<float>();
  mytree->Branch ("muonSize", &m_muonSize);
  m_Z_mass = new std::vector<float>();
  mytree->Branch ("Z_mass", &m_Z_mass);
  m_truthZ = new std::vector<float>();
  mytree->Branch ("truthZ", &m_truthZ);
  m_children = new std::vector<float>();
  mytree->Branch ("children", &m_children);
  m_mChildren = new std::vector<float>();
  mytree->Branch ("mChildren", &m_mChildren);
  m_Z_match = new std::vector<float>();
  mytree->Branch ("Z_match", &m_Z_match);
  m_truthEta = new std::vector<float>();
  mytree->Branch ("truthEta", &m_truthEta);
  m_truthPhi = new std::vector<float>();
  mytree->Branch ("truthPhi", &m_truthPhi);
  m_matchEta = new std::vector<float>();
  mytree->Branch ("matchEta", &m_matchEta);
  m_matchPhi = new std::vector<float>();
  mytree->Branch ("matchPhi", &m_matchPhi);
  m_ptRes = new std::vector<float>();
  mytree->Branch ("ptRes", &m_ptRes);
  m_etaRes = new std::vector<float>();
  mytree->Branch ("etaRes", &m_etaRes);
  m_MinBias = new std::vector<float>();
  mytree->Branch ("MinBias", &m_MinBias);
  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  ANA_MSG_INFO ("in execute");
  m_muonSize->clear();
  m_muonEta->clear();
  m_muonPhi->clear();
  m_muonPt->clear();
  m_muonE->clear();
  m_leadingPt->clear();
  m_subPt->clear();
  m_Z_mass->clear();
  m_truthZ->clear();
  m_children->clear();
  m_truthEta->clear();
  m_truthPhi->clear();
  m_mChildren->clear();
  m_Z_match->clear();
  m_MinBias->clear();
  m_matchEta->clear();
  m_matchPhi->clear();
  m_ptRes->clear();
  m_etaRes->clear();
  const xAOD::EventInfo* eventInfo = 0;
  ANA_CHECK(evtStore()->retrieve( eventInfo, "EventInfo"));
  // check if the event is data or MC
  // (many tools are applied either to data or MC)
  bool isMC = false;
  // check if the event is MC
  if (eventInfo->eventType (xAOD::EventInfo::IS_SIMULATION)) {
    isMC = true; // can do something with this later
    ANA_MSG_INFO ("MC is triggering");
    const xAOD::TruthEventContainer* xTruthEventContainer = NULL;
    ANA_CHECK(evtStore()->retrieve(xTruthEventContainer, "TruthEvents"));
    xAOD::TruthEventContainer::const_iterator itr;
    auto tHSevent = xTruthEventContainer->at(0);
    int nPart = tHSevent->nTruthParticles();
    for (int iPart = 0; iPart < nPart; iPart++) {
      const xAOD::TruthParticle* particle = tHSevent->truthParticle(iPart);
      if (particle) {
        if (particle->pdgId() == 23) {
          double_t truthZ = (particle->p4()).M();
          const xAOD::TruthParticle* child0 = particle->child(0);
          const xAOD::TruthParticle* child1 = particle->child(1);
          if (child0) {
            if (child1) {
              if (child0->pdgId() == 13 && child1->pdgId() == -13) {
                m_truthZ->push_back (truthZ * 0.001);
                m_children->push_back (child0->pt() * 0.001);
                m_children->push_back (child1->pt() * 0.001);
                m_truthEta->push_back (child0->eta());
                m_truthEta->push_back (child1->eta());
                m_truthPhi->push_back (child0->phi());
                m_truthPhi->push_back (child1->phi()); }}}}}}
    const xAOD::MuonContainer* muons = nullptr;
    ANA_CHECK (evtStore()->retrieve (muons, "Muons"));
    m_muonSize->push_back (muons->size());
    if (muons->size() < 2) {          //ignore events with less than 2 muons
      ANA_MSG_INFO ("No Muons Here");
      tree ("analysis")->Fill ();
      return StatusCode::SUCCESS; }
    const xAOD::Muon* muon_m1 = nullptr;
    const xAOD::Muon* muon_m2 = nullptr;
    const xAOD::Muon* MinBias1 = nullptr;
    const xAOD::Muon* MinBias2 = nullptr;
    Double_t ptRes1 = 0;
    Double_t ptRes2 = 0;
    Double_t etaRes1 = 0;
    Double_t etaRes2 = 0;
    Double_t value1 = 0;
    Double_t value2 = 0;
    for (const xAOD::Muon* muon : *muons) {
      typedef ElementLink<xAOD::TruthParticleContainer>ElementTruthLink_t;
      const xAOD::TruthParticle* tresult = 0;
      if (muon->isAvailable<ElementTruthLink_t>("truthParticleLink")) {
        const ElementTruthLink_t ptruthContainer = muon->auxdata<ElementTruthLink_t>("truthParticleLink");
        if (ptruthContainer.isValid()) {
          tresult = *ptruthContainer; }}
      if (tresult != 0) {
        int pdgId = tresult->pdgId();
        if (pdgId == 13 || pdgId == -13) {
          const xAOD::TruthParticle* mother = tresult->parent(0);
          if (mother) {
          int pdgIdParent = mother->pdgId();
            if (pdgIdParent == 23) {
              if (muon_m1 == nullptr) {
                muon_m1 = muon;
                ptRes1 = ((1/(tresult->pt() * 0.001)) - (1/(muon_m1->pt() * 0.001)));
                etaRes1 = (tresult->eta() - muon_m1->eta()); }
                else {
                muon_m2 = muon;
                ptRes2 = ((1/(tresult->pt() * 0.001)) - (1/(muon_m2->pt() * 0.001)));
                etaRes2 = (tresult->eta() - muon_m2->eta()); }}}}}
      if (tresult == 0) {
        if (muon->pt() > value1) {
          value1 = muon->pt();
          MinBias1 = muon; }}}
    for (const xAOD::Muon* muon : *muons) {
      typedef ElementLink<xAOD::TruthParticleContainer>ElementTruthLink_t;
      const xAOD::TruthParticle* tresult = 0;
      if (muon->isAvailable<ElementTruthLink_t>("truthParticleLink")) {
        const ElementTruthLink_t ptruthContainer = muon->auxdata<ElementTruthLink_t>("truthParticleLink");
        if (ptruthContainer.isValid()) {
          tresult = *ptruthContainer; }}
      if (tresult == 0) {
        if (muon->pt() > value2) {
          if (muon->pt() == value1) {continue;}
            value2 = muon->pt();
            MinBias2 = muon; }}}
    if (muon_m1 != nullptr) {
      m_mChildren->push_back (muon_m1->pt() * 0.001);
      m_matchEta->push_back (muon_m1->eta());
      m_matchPhi->push_back (muon_m1->phi());
      m_ptRes->push_back (ptRes1);
      m_etaRes->push_back (etaRes1);
      if (muon_m2 != nullptr) {
        m_mChildren->push_back (muon_m2->pt() * 0.001);
        m_matchEta->push_back (muon_m2->eta());
        m_matchPhi->push_back (muon_m2->phi());
        m_ptRes->push_back (ptRes2);
        m_etaRes->push_back (etaRes2);
        TLorentzVector muon_m1_p4 = muon_m1->p4();
        TLorentzVector muon_m2_p4 = muon_m2->p4();
        TLorentzVector Zm_p4 = muon_m1_p4 + muon_m2_p4;
        Double_t Z_match = Zm_p4.M();
        m_Z_match->push_back (Z_match * 0.001); }}
    if (MinBias1 != nullptr) {
      if (MinBias2 != nullptr) {
        TLorentzVector MinBias1_p4 = MinBias1->p4();
        TLorentzVector MinBias2_p4 = MinBias2->p4();
        TLorentzVector MB_p4 = MinBias1_p4 + MinBias2_p4;
        Double_t MinBias = MB_p4.M();
        m_MinBias->push_back (MinBias * 0.001); }}
    Double_t pt_value = 0;
    Double_t pt_sub = 0;
    const xAOD::Muon* muon_leading = nullptr;
    const xAOD::Muon* muon_sub = nullptr;
    for (const xAOD::Muon* muon : *muons) {
      ANA_MSG_INFO ("execute(): original muon pt = " << ((muon)->pt() * 0.001) << " GeV");
      if (muon->pt() > pt_value) {                    //find leading muon
        pt_value = muon->pt();
        muon_leading = muon; } }
    for (const xAOD::Muon* muon : *muons) {
      if (muon->pt() > pt_sub) {                      //find sub-leading muon
        if (muon->pt() == pt_value) {continue;}       //ignore leading muon
          pt_sub = muon->pt();
          muon_sub = muon; }}
    TLorentzVector muon_leading_p4 = muon_leading->p4(); //define muon 4-momenta
    TLorentzVector muon_sub_p4 = muon_sub->p4();
    TLorentzVector Z_p4 = muon_leading_p4 + muon_sub_p4;
    Double_t Z_mass = Z_p4.M();                          //calculate Z mass
    const xAOD::EventInfo* ei = nullptr;
    ANA_CHECK (evtStore()->retrieve (ei, "EventInfo"));
    m_runNumber = ei->runNumber ();
    m_eventNumber = ei->eventNumber ();
    for (const xAOD::Muon* muon : *muons) {
      m_muonEta->push_back (muon->eta());
      m_muonPhi->push_back (muon->phi());
      m_muonPt->push_back (muon->pt() * 0.001);   //factor 0.001 to make unit GeV
      m_muonE->push_back (muon->e() * 0.001); }
    m_leadingPt->push_back (muon_leading->pt() * 0.001);
    m_subPt->push_back (muon_sub->pt() * 0.001);
    m_Z_mass->push_back (Z_mass * 0.001);
    tree ("analysis")->Fill ();
  }
  // if data check if event passes GRL
  if (isMC == false) { // it's data!
    ANA_MSG_INFO ("data should trigger");
    if (!m_grl->passRunLB(*eventInfo)) {
      ANA_MSG_INFO ("drop event: GRL");
      return StatusCode::SUCCESS; } // go to next event
    const xAOD::MuonContainer* muons = nullptr;
    ANA_CHECK (evtStore()->retrieve (muons, "AnalysisMuons_NOSYS"));
    m_muonSize->push_back (muons->size());
    if (muons->size() < 2) {          //ignore events with less than 2 muons
      ANA_MSG_INFO ("No Muons Here");
      tree ("analysis")->Fill ();
      return StatusCode::SUCCESS; }
    Double_t pt_value = 0;
    Double_t pt_sub = 0;
    const xAOD::Muon* muon_leading = nullptr;
    const xAOD::Muon* muon_sub = nullptr;
    for (const xAOD::Muon* muon : *muons) {
      ANA_MSG_INFO ("execute(): original muon pt = " << ((muon)->pt() * 0.001) << " GeV");
      if (muon->pt() > pt_value) {                    //find leading muon
        pt_value = muon->pt();
        muon_leading = muon; } }
    for (const xAOD::Muon* muon : *muons) {
      if (muon->pt() > pt_sub) {                      //find sub-leading muon
        if (muon->pt() == pt_value) {continue;}       //ignore leading muon
          pt_sub = muon->pt();
          muon_sub = muon; } 
    }
    TLorentzVector muon_leading_p4 = muon_leading->p4(); //define muon 4-momenta
    TLorentzVector muon_sub_p4 = muon_sub->p4();
    TLorentzVector Z_p4 = muon_leading_p4 + muon_sub_p4;
    Double_t Z_mass = Z_p4.M();                          //calculate Z mass
    const xAOD::EventInfo* ei = nullptr;
    ANA_CHECK (evtStore()->retrieve (ei, "EventInfo_NOSYS"));
    m_runNumber = ei->runNumber ();
    m_eventNumber = ei->eventNumber ();
    for (const xAOD::Muon* muon : *muons) {
      m_muonEta->push_back (muon->eta());
      m_muonPhi->push_back (muon->phi());
      m_muonPt->push_back (muon->pt() * 0.001);   //factor 0.001 to make unit GeV
      m_muonE->push_back (muon->e() * 0.001); }
    m_leadingPt->push_back (muon_leading->pt() * 0.001);
    m_subPt->push_back (muon_sub->pt() * 0.001);
    m_Z_mass->push_back (Z_mass * 0.001);
    tree ("analysis")->Fill ();
  } // end if not MC
  ANA_MSG_INFO ("keep event: GRL");
  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.
  return StatusCode::SUCCESS;
}
