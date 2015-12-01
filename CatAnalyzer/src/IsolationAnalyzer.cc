#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include <vector>

using namespace pat;
using namespace edm;

class IsolationAnalyzer : public edm::EDAnalyzer {
public:
  explicit IsolationAnalyzer(const edm::ParameterSet&);
  ~IsolationAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexLabel_;
  edm::EDGetTokenT<edm::View<pat::Muon> >          muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> >      electronToken_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate> >  packedCandidateToken_;
  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

  //double GetIsolation( const reco::Candidate::LorentzVector& pasObj, edm::Handle<edm::View<pat::PackedCandidate> > packedCandidates, const double dR);
  std::vector<double> GetIsolation( const reco::Candidate::LorentzVector& pasObj, edm::Handle<edm::View<pat::PackedCandidate> > packedCandidates, const double dR);
  bool isGenIdMatched( const reco::Candidate::LorentzVector& lepton, edm::Handle<reco::GenParticleCollection> genParticles_, int id );
  bool isFromWorZ( const reco::Candidate::LorentzVector& lepton, Handle<reco::GenParticleCollection> genParticles_, int id );
  bool MatchObjects( const reco::Candidate::LorentzVector& pasObj,
      const reco::Candidate::LorentzVector& proObj,
      bool exact );

  TTree *tree;
  float weight;
  float mu1_pt;
  float mu1_eta;
  float mu1_pfRelIso;
  float mu1_chIso;
  float mu1_nhIso;
  float mu1_phIso;
  float mu1_puChIso;

  float mu1_pfRelIso03;
  float mu1_chIso03;
  float mu1_nhIso03;
  float mu1_phIso03;
  float mu1_puChIso03;

  float mu1_pfRelIso04;
  float mu1_chIso04;
  float mu1_nhIso04;
  float mu1_phIso04;
  float mu1_puChIso04;

  float mu1_pfRelIso03_packed;
  float mu1_chIso03_packed;
  float mu1_nhIso03_packed;
  float mu1_phIso03_packed;
  float mu1_puChIso03_packed;

  float mu1_pfRelIso04_packed;
  float mu1_chIso04_packed;
  float mu1_nhIso04_packed;
  float mu1_phIso04_packed;
  float mu1_puChIso04_packed;

  float mu1_detRelIso;
  float mu1_trackIso;
  float mu1_ecalIso;
  float mu1_hcalIso;
  int mu1_tightID;
  int mu1_matched;
  int mu1_matchedZ;

  float mu2_pt;
  float mu2_eta;
  float mu2_pfRelIso;
  float mu2_chIso;
  float mu2_nhIso;
  float mu2_phIso;
  float mu2_puChIso;

  float mu2_pfRelIso03;
  float mu2_chIso03;
  float mu2_nhIso03;
  float mu2_phIso03;
  float mu2_puChIso03;

  float mu2_pfRelIso04;
  float mu2_chIso04;
  float mu2_nhIso04;
  float mu2_phIso04;
  float mu2_puChIso04;

  float mu2_pfRelIso03_packed;
  float mu2_chIso03_packed;
  float mu2_nhIso03_packed;
  float mu2_phIso03_packed;
  float mu2_puChIso03_packed;

  float mu2_pfRelIso04_packed;
  float mu2_chIso04_packed;
  float mu2_nhIso04_packed;
  float mu2_phIso04_packed;
  float mu2_puChIso04_packed;

  float mu2_detRelIso;
  float mu2_trackIso;
  float mu2_ecalIso;
  float mu2_hcalIso;
  int mu2_tightID;
  int mu2_matched;
  int mu2_matchedZ;


  float el1_pt;
  float el1_eta;
  float el1_chIso;
  float el1_nhIso;
  float el1_phIso;
  float el1_pfRelIso;
  float el1_detRelIso;
  float el1_trackIso;
  float el1_ecalIso;
  float el1_hcalIso;
  int el1_matched;
  int el1_matchedZ;
  int el1_mediumID;
  int el1_tightID;
 
};

IsolationAnalyzer::IsolationAnalyzer(const edm::ParameterSet& iConfig)
{
  genToken_      = consumes<reco::GenParticleCollection>   (iConfig.getParameter<edm::InputTag>("genLabel"));
  vertexLabel_      = consumes<reco::VertexCollection>       (iConfig.getParameter<edm::InputTag>("vertexLabel"));
  muonToken_     = consumes<edm::View<pat::Muon> >          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_ = consumes<edm::View<pat::Electron> >      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  packedCandidateToken_ = consumes<edm::View<pat::PackedCandidate> >      (iConfig.getParameter<edm::InputTag>("packedCandiateLabel"));
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >     (iConfig.getParameter<edm::InputTag>("eleMediumIdMap"));
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >      (iConfig.getParameter<edm::InputTag>("eleTightIdMap"));

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("weight", &weight, "weight/F");
  tree->Branch("mu1_pt", &mu1_pt, "mu1_pt/F");
  tree->Branch("mu1_eta", &mu1_eta, "mu1_eta/F");
  tree->Branch("mu1_pfRelIso", &mu1_pfRelIso, "mu1_pfRelIso/F");
  tree->Branch("mu1_chIso", &mu1_chIso, "mu1_chIso/F");
  tree->Branch("mu1_nhIso", &mu1_nhIso, "mu1_nhIso/F");
  tree->Branch("mu1_phIso", &mu1_phIso, "mu1_phIso/F");
  tree->Branch("mu1_puChIso", &mu1_puChIso, "mu1_puChIso/F");

  tree->Branch("mu1_pfRelIso03", &mu1_pfRelIso03, "mu1_pfRelIso03/F");
  tree->Branch("mu1_chIso03", &mu1_chIso03, "mu1_chIso03/F");
  tree->Branch("mu1_nhIso03", &mu1_nhIso03, "mu1_nhIso03/F");
  tree->Branch("mu1_phIso03", &mu1_phIso03, "mu1_phIso03/F");
  tree->Branch("mu1_puChIso03", &mu1_phIso03, "mu1_phIso03/F");

  tree->Branch("mu1_pfRelIso04", &mu1_pfRelIso04, "mu1_pfRelIso04/F");
  tree->Branch("mu1_chIso04", &mu1_chIso04, "mu1_chIso04/F");
  tree->Branch("mu1_nhIso04", &mu1_nhIso04, "mu1_nhIso04/F");
  tree->Branch("mu1_phIso04", &mu1_phIso04, "mu1_phIso04/F");
  tree->Branch("mu1_puChIso04", &mu1_phIso04, "mu1_phIso04/F");
 
  tree->Branch("mu1_pfRelIso03_packed", &mu1_pfRelIso03_packed, "mu1_pfRelIso03_packed/F");
  tree->Branch("mu1_chIso03_packed", &mu1_chIso03_packed, "mu1_chIso03_packed/F");
  tree->Branch("mu1_nhIso03_packed", &mu1_nhIso03_packed, "mu1_nhIso03_packed/F");
  tree->Branch("mu1_phIso03_packed", &mu1_phIso03_packed, "mu1_phIso03_packed/F");
  tree->Branch("mu1_puChIso03_packed", &mu1_phIso03_packed, "mu1_phIso03_packed/F");

  tree->Branch("mu1_pfRelIso04_packed", &mu1_pfRelIso04_packed, "mu1_pfRelIso04_packed/F");
  tree->Branch("mu1_chIso04_packed", &mu1_chIso04_packed, "mu1_chIso04_packed/F");
  tree->Branch("mu1_nhIso04_packed", &mu1_nhIso04_packed, "mu1_nhIso04_packed/F");
  tree->Branch("mu1_phIso04_packed", &mu1_phIso04_packed, "mu1_phIso04_packed/F");
  tree->Branch("mu1_puChIso04_packed", &mu1_phIso04_packed, "mu1_phIso04_packed/F");

  tree->Branch("mu1_detRelIso", &mu1_detRelIso, "mu1_detRelIso/F");
  tree->Branch("mu1_trackIso", &mu1_trackIso, "mu1_trackIso/F");
  tree->Branch("mu1_ecalIso", &mu1_ecalIso, "mu1_ecalIso/F");
  tree->Branch("mu1_hcalIso", &mu1_hcalIso, "mu1_hcalIso/F");
  tree->Branch("mu1_tightID", &mu1_tightID, "mu1_tightID/I");
  tree->Branch("mu1_matched", &mu1_matched, "mu1_matched/I");
  tree->Branch("mu1_matchedZ", &mu1_matchedZ, "mu1_matchedZ/I");

  tree->Branch("mu2_pt", &mu2_pt, "mu2_pt/F");
  tree->Branch("mu2_eta", &mu2_eta, "mu2_eta/F");
  tree->Branch("mu2_pfRelIso", &mu2_pfRelIso, "mu2_pfRelIso/F");
  tree->Branch("mu2_chIso", &mu2_chIso, "mu2_chIso/F");
  tree->Branch("mu2_nhIso", &mu2_nhIso, "mu2_nhIso/F");
  tree->Branch("mu2_phIso", &mu2_phIso, "mu2_phIso/F");
  tree->Branch("mu2_puChIso", &mu2_puChIso, "mu2_puChIso/F");

  tree->Branch("mu2_pfRelIso03", &mu2_pfRelIso03, "mu2_pfRelIso03/F");
  tree->Branch("mu2_chIso03", &mu2_chIso03, "mu2_chIso03/F");
  tree->Branch("mu2_nhIso03", &mu2_nhIso03, "mu2_nhIso03/F");
  tree->Branch("mu2_phIso03", &mu2_phIso03, "mu2_phIso03/F");
  tree->Branch("mu2_puChIso03", &mu2_phIso03, "mu2_phIso03/F");

  tree->Branch("mu2_pfRelIso04", &mu2_pfRelIso04, "mu2_pfRelIso04/F");
  tree->Branch("mu2_chIso04", &mu2_chIso04, "mu2_chIso04/F");
  tree->Branch("mu2_nhIso04", &mu2_nhIso04, "mu2_nhIso04/F");
  tree->Branch("mu2_phIso04", &mu2_phIso04, "mu2_phIso04/F");
  tree->Branch("mu2_puChIso04", &mu2_phIso04, "mu2_phIso04/F");

  tree->Branch("mu2_pfRelIso03_packed", &mu2_pfRelIso03_packed, "mu2_pfRelIso03_packed/F");
  tree->Branch("mu2_chIso03_packed", &mu2_chIso03_packed, "mu2_chIso03_packed/F");
  tree->Branch("mu2_nhIso03_packed", &mu2_nhIso03_packed, "mu2_nhIso03_packed/F");
  tree->Branch("mu2_phIso03_packed", &mu2_phIso03_packed, "mu2_phIso03_packed/F");
  tree->Branch("mu2_puChIso03_packed", &mu2_phIso03_packed, "mu2_phIso03_packed/F");

  tree->Branch("mu2_pfRelIso04_packed", &mu2_pfRelIso04_packed, "mu2_pfRelIso04_packed/F");
  tree->Branch("mu2_chIso04_packed", &mu2_chIso04_packed, "mu2_chIso04_packed/F");
  tree->Branch("mu2_nhIso04_packed", &mu2_nhIso04_packed, "mu2_nhIso04_packed/F");
  tree->Branch("mu2_phIso04_packed", &mu2_phIso04_packed, "mu2_phIso04_packed/F");
  tree->Branch("mu2_puChIso04_packed", &mu2_phIso04_packed, "mu2_phIso04_packed/F");

  tree->Branch("mu2_detRelIso", &mu2_detRelIso, "mu2_detRelIso/F");
  tree->Branch("mu2_trackIso", &mu2_trackIso, "mu2_trackIso/F");
  tree->Branch("mu2_ecalIso", &mu2_ecalIso, "mu2_ecalIso/F");
  tree->Branch("mu2_hcalIso", &mu2_hcalIso, "mu2_hcalIso/F");
  tree->Branch("mu2_tightID", &mu2_tightID, "mu2_tightID/I");
  tree->Branch("mu2_matched", &mu2_matched, "mu2_matched/I");
  tree->Branch("mu2_matchedZ", &mu2_matchedZ, "mu2_matchedZ/I");

  tree->Branch("el1_pt", &el1_pt, "el1_pt/F");
  tree->Branch("el1_eta", &el1_eta, "el1_eta/F");
  tree->Branch("el1_pfRelIso", &el1_pfRelIso, "el1_pfRelIso/F");
  tree->Branch("el1_chIso", &el1_chIso, "el1_chIso/F");
  tree->Branch("el1_nhIso", &el1_nhIso, "el1_nhIso/F");
  tree->Branch("el1_phIso", &el1_phIso, "el1_phIso/F");
  tree->Branch("el1_detRelIso", &el1_detRelIso, "el1_detRelIso/F");
  tree->Branch("el1_trackIso", &el1_pfRelIso, "el1_trackIso/F");
  tree->Branch("el1_ecalIso", &el1_ecalIso, "el1_ecalIso/F");
  tree->Branch("el1_hcalIso", &el1_hcalIso, "el1_hcalIso/F");
  tree->Branch("el1_matched", &el1_matched, "el1_matched/I");
  tree->Branch("el1_matchedZ", &el1_matchedZ, "el1_matchedZ/I");
  tree->Branch("el1_mediumID", &el1_mediumID, "el1_mediumID/I");
  tree->Branch("el1_tightID", &el1_mediumID, "el1_tightID/I");

}


IsolationAnalyzer::~IsolationAnalyzer()
{

}

void IsolationAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  weight = 1.0;

  mu1_pt = -1.0;
  mu1_eta = -9.0;
  mu1_pfRelIso = -1.0;
  mu1_chIso = -1.0;
  mu1_nhIso = -1.0;
  mu1_phIso = -1.0;
  mu1_puChIso = -1.0;

  mu1_pfRelIso03 = -1.0;
  mu1_chIso03 = -1.0;
  mu1_nhIso03 = -1.0;
  mu1_phIso03 = -1.0;
  mu1_puChIso03 = -1.0;

  mu1_pfRelIso04 = -1.0;
  mu1_chIso04 = -1.0;
  mu1_nhIso04 = -1.0;
  mu1_phIso04 = -1.0;
  mu1_puChIso04 = -1.0;

  mu1_pfRelIso03_packed = -1.0;
  mu1_chIso03_packed = -1.0;
  mu1_nhIso03_packed = -1.0;
  mu1_phIso03_packed = -1.0;
  mu1_puChIso03_packed = -1.0;

  mu1_pfRelIso04_packed = -1.0;
  mu1_chIso04_packed = -1.0;
  mu1_nhIso04_packed = -1.0;
  mu1_phIso04_packed = -1.0;
  mu1_puChIso04_packed = -1.0;

  mu1_detRelIso = -1.0;
  mu1_trackIso = -1.0;
  mu1_hcalIso = -1.0;
  mu1_ecalIso = -1.0;
  mu1_tightID = -1.0;
  mu1_matched = -1.0;
  mu1_matchedZ = -1.0;

  mu2_pt = -1.0;
  mu2_eta = -9.0;
  mu2_pfRelIso = -1.0;
  mu2_chIso = -1.0;
  mu2_nhIso = -1.0;
  mu2_phIso = -1.0;
  mu2_puChIso = -1.0;

  mu2_pfRelIso03 = -1.0;
  mu2_chIso03 = -1.0;
  mu2_nhIso03 = -1.0;
  mu2_phIso03 = -1.0;
  mu2_puChIso03 = -1.0;

  mu2_pfRelIso04 = -1.0;
  mu2_chIso04 = -1.0;
  mu2_nhIso04 = -1.0;
  mu2_phIso04 = -1.0;
  mu2_puChIso04 = -1.0;

  mu2_pfRelIso03_packed = -1.0;
  mu2_chIso03_packed = -1.0;
  mu2_nhIso03_packed = -1.0;
  mu2_phIso03_packed = -1.0;
  mu2_puChIso03_packed = -1.0;

  mu2_pfRelIso04_packed = -1.0;
  mu2_chIso04_packed = -1.0;
  mu2_nhIso04_packed = -1.0;
  mu2_phIso04_packed = -1.0;
  mu2_puChIso04_packed = -1.0;

  mu2_detRelIso = -1.0;
  mu2_trackIso = -1.0;
  mu2_hcalIso = -1.0;
  mu2_ecalIso = -1.0;
  mu2_tightID = -1.0;
  mu2_matched = -1.0;
  mu2_matchedZ = -1.0;

  el1_pt = -1.0;
  el1_eta = -1.0;
  el1_pfRelIso = -1.0;
  el1_chIso = -1.0;
  el1_nhIso = -1.0;
  el1_phIso = -1.0;
  el1_trackIso = -1.0;
  el1_hcalIso = -1.0;
  el1_ecalIso = -1.0;
  el1_matched = -1;
  el1_matchedZ = -1;
  el1_mediumID = -1;
  el1_tightID = -1;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genToken_, genParticles);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByToken(vertexLabel_,recVtxs);
  reco::Vertex pv;
  if (recVtxs->size())
    pv = recVtxs->at(0);


  Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(electronToken_, electrons);

  Handle<edm::View<pat::PackedCandidate> > packedCandidates;
  iEvent.getByToken(packedCandidateToken_, packedCandidates);

   // Get the electron ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);

  for (unsigned int i = 0; i < electrons->size() ; i++) {
    const pat::Electron & electron = electrons->at(i);

    bool id = electron.pt() > 20 && fabs( electron.eta() ) < 2.4 && electron.isPF() ;

    const auto el = electrons->ptrAt(i);
    bool isPassMedium = (*medium_id_decisions)[el];
    bool isPassTight  = (*tight_id_decisions)[el];

    if( id == false ) continue;

    el1_pt = electron.pt();
    el1_eta = electron.eta();
    el1_chIso = electron.chargedHadronIso();
    el1_nhIso = electron.neutralHadronIso();
    el1_phIso = electron.photonIso();
    el1_pfRelIso = ( el1_chIso + el1_nhIso + el1_phIso ) / el1_pt ;
    el1_trackIso = electron.trackIso();
    el1_ecalIso = electron.ecalIso();
    el1_hcalIso = electron.hcalIso();
    el1_detRelIso = ( el1_trackIso + el1_ecalIso + el1_hcalIso ) / el1_pt;
    el1_matched = isGenIdMatched( electron.p4(), genParticles, 11);
    el1_matchedZ = isFromWorZ( electron.p4(), genParticles, 11);
    el1_mediumID = isPassMedium;
    el1_tightID = isPassTight;
  }

  Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonToken_, muons);

  int idx = 0;
  for (unsigned int i = 0; i < muons->size() ; i++) {
    const pat::Muon & muon = muons->at(i);
    bool id = muon.pt() > 20 && fabs( muon.eta() ) < 2.4 && muon.isPFMuon() ;
    bool tightid = muon.isTightMuon(pv);

    if( id == false ) continue;
    if( idx == 0){
      mu1_pt = muon.pt();
      mu1_eta= muon.eta();
      mu1_chIso = muon.chargedHadronIso();
      mu1_nhIso = muon.neutralHadronIso();
      mu1_phIso = muon.photonIso();
      mu1_puChIso = muon.puChargedHadronIso();
      mu1_pfRelIso = ( mu1_chIso + mu1_nhIso + mu1_phIso ) / mu1_pt ;

      mu1_chIso04 = muon.pfIsolationR04().sumChargedHadronPt;
      mu1_nhIso04 = muon.pfIsolationR04().sumNeutralHadronEt;
      mu1_phIso04 = muon.pfIsolationR04().sumPhotonEt;
      mu1_chIso03 = muon.pfIsolationR03().sumChargedHadronPt;
      mu1_nhIso03 = muon.pfIsolationR03().sumNeutralHadronEt;
      mu1_phIso03 = muon.pfIsolationR03().sumPhotonEt;

      mu1_pfRelIso03 = ( mu1_chIso03 + mu1_nhIso03 + mu1_phIso03 ) / mu1_pt ;
      mu1_pfRelIso04 = ( mu1_chIso04 + mu1_nhIso04 + mu1_phIso04 ) / mu1_pt ;

      std::vector<double> Iso03_packedCandidate = GetIsolation( muon.p4(), packedCandidates, 0.3);
      std::vector<double> Iso04_packedCandidate = GetIsolation( muon.p4(), packedCandidates, 0.4);

      mu1_chIso03_packed = Iso03_packedCandidate[0];
      mu1_chIso04_packed = Iso04_packedCandidate[0];
 
      mu1_nhIso03_packed = Iso03_packedCandidate[1];
      mu1_nhIso04_packed = Iso04_packedCandidate[1];

      mu1_phIso03_packed = Iso03_packedCandidate[2];
      mu1_phIso04_packed = Iso04_packedCandidate[2];

      mu1_pfRelIso03_packed =  ( Iso03_packedCandidate[0] + Iso03_packedCandidate[1] + Iso03_packedCandidate[2]) / mu1_pt ;
      mu1_pfRelIso04_packed =  ( Iso04_packedCandidate[0] + Iso04_packedCandidate[1] + Iso04_packedCandidate[2]) / mu1_pt ;

      //std::cout << "packed   isolation 0.3 = " << Iso03_packedCandidate[0] << " " << Iso03_packedCandidate[1] << " " << Iso03_packedCandidate[2]  << std::endl;
      //std::cout << "function isolation 0.3 = " << mu1_chIso03 << " " << mu1_nhIso03 << " " << mu1_phIso03  << std::endl;

      mu1_trackIso = muon.trackIso();
      mu1_ecalIso = muon.ecalIso();
      mu1_hcalIso = muon.hcalIso();
      mu1_detRelIso = ( mu1_trackIso + mu1_ecalIso + mu1_hcalIso ) / mu1_pt;
      mu1_tightID = tightid;
      mu1_matched = isGenIdMatched( muon.p4(), genParticles, 13);
      mu1_matchedZ = isFromWorZ( muon.p4(), genParticles, 13);
    }
    if( idx == 1){
      mu2_pt = muon.pt();
      mu2_eta = muon.eta();
      mu2_chIso = muon.chargedHadronIso();
      mu2_nhIso = muon.neutralHadronIso();
      mu2_phIso = muon.photonIso();
      mu2_puChIso = muon.puChargedHadronIso();      
      mu2_pfRelIso = ( mu2_chIso + mu2_nhIso + mu2_phIso ) / mu2_pt ;

      mu2_chIso04 = muon.pfIsolationR04().sumChargedHadronPt;
      mu2_nhIso04 = muon.pfIsolationR04().sumNeutralHadronEt;
      mu2_phIso04 = muon.pfIsolationR04().sumPhotonEt;
      mu2_chIso03 = muon.pfIsolationR03().sumChargedHadronPt;
      mu2_nhIso03 = muon.pfIsolationR03().sumNeutralHadronEt;
      mu2_phIso03 = muon.pfIsolationR03().sumPhotonEt;

      mu2_pfRelIso03 = ( mu2_chIso03 + mu2_nhIso03 + mu2_phIso03 ) / mu2_pt ;
      mu2_pfRelIso04 = ( mu2_chIso04 + mu2_nhIso04 + mu2_phIso04 ) / mu2_pt ;

      std::vector<double> Iso03_packedCandidate = GetIsolation( muon.p4(), packedCandidates, 0.3);
      std::vector<double> Iso04_packedCandidate = GetIsolation( muon.p4(), packedCandidates, 0.4);

      mu2_chIso03_packed = Iso03_packedCandidate[0];
      mu2_chIso04_packed = Iso04_packedCandidate[0];
 
      mu2_nhIso03_packed = Iso03_packedCandidate[1];
      mu2_nhIso04_packed = Iso04_packedCandidate[1];

      mu2_phIso03_packed = Iso03_packedCandidate[2];
      mu2_phIso04_packed = Iso04_packedCandidate[2];

      mu2_pfRelIso03_packed =  ( Iso03_packedCandidate[0] + Iso03_packedCandidate[1] + Iso03_packedCandidate[2]) / mu2_pt ;
      mu2_pfRelIso04_packed =  ( Iso04_packedCandidate[0] + Iso04_packedCandidate[1] + Iso04_packedCandidate[2]) / mu2_pt ;

      mu2_trackIso = muon.trackIso();
      mu2_ecalIso = muon.ecalIso();
      mu2_hcalIso = muon.hcalIso();
      mu2_detRelIso = ( mu2_trackIso + mu2_ecalIso + mu2_hcalIso ) / mu2_pt;
      mu2_tightID = tightid;
      mu2_matched = isGenIdMatched( muon.p4(), genParticles, 13);
      mu2_matchedZ = isFromWorZ( muon.p4(), genParticles, 13);
    }
    idx++;
   
  }

  tree->Fill();
}

std::vector<double> IsolationAnalyzer::GetIsolation( const reco::Candidate::LorentzVector& pasObj, edm::Handle<edm::View<pat::PackedCandidate> > packedCandidates, const double dR ){

  double chIso = 0;
  double nhIso = 0;
  double phIso = 0;

  for( edm::View<pat::PackedCandidate>::const_iterator mcIter=packedCandidates->begin(); mcIter != packedCandidates->end(); mcIter++){
    
    double proEta = mcIter->eta();
    double proPhi = mcIter->phi();
    double proPt  = mcIter->pt();
    double pasEta = pasObj.eta();
    double pasPhi = pasObj.phi();
    //double pasPt  = pasObj.pt();

    double dRval = deltaR(proEta, proPhi, pasEta, pasPhi);

    if( abs(mcIter->pdgId()) == 211 && dRval < dR){
      chIso = chIso + proPt;
    } 

    if( abs(mcIter->pdgId()) == 130 && dRval < dR){
      nhIso = chIso + proPt;
    }

    if( abs(mcIter->pdgId()) == 22 && dRval < dR){
      phIso = chIso + proPt;
    }

  }

  std::vector<double>  output;
  output.push_back(chIso);
  output.push_back(nhIso);
  output.push_back(phIso);
  //double absIso = chIso + nhIso + phIso;
  //double relIso = absIso / pasObj.pt();
  //return relIso;
  return output;

}

bool IsolationAnalyzer::isGenIdMatched( const reco::Candidate::LorentzVector& lepton, edm::Handle<reco::GenParticleCollection> genParticles_, int id ){

  bool out = false;

  for ( reco::GenParticleCollection::const_iterator mcIter=genParticles_->begin(); mcIter != genParticles_->end(); mcIter++ ) {
    int genId = mcIter->pdgId();
    if( abs(genId) != (int) id ) continue;
      
    bool match = MatchObjects(lepton, mcIter->p4(), false); 
    if( match != true) continue;

    out = true;
  }

  return out;
}

bool IsolationAnalyzer::isFromWorZ( const reco::Candidate::LorentzVector& lepton, edm::Handle<reco::GenParticleCollection> genParticles_, int id ){

  bool out = false;

  for ( reco::GenParticleCollection::const_iterator mcIter=genParticles_->begin(); mcIter != genParticles_->end(); mcIter++ ) {
    int genId = mcIter->pdgId();

    if( abs(genId) != (int) id ) continue;

    bool match = MatchObjects(lepton, mcIter->p4(), false);

    if( match != true) continue;
   
    const reco::Candidate* mother = mcIter->mother();
    while( mother != 0 ){
      if( abs(mother->pdgId()) == 24 || abs(mother->pdgId()) == 23 ) { 
        out = true;
      }
      mother = mother->mother();
    }
  }

  return out;
}

bool IsolationAnalyzer::MatchObjects( const reco::Candidate::LorentzVector& pasObj,
      const reco::Candidate::LorentzVector& proObj,
      bool exact ) {
    double proEta = proObj.eta();
    double proPhi = proObj.phi();
    double proPt  = proObj.pt();
    double pasEta = pasObj.eta();
    double pasPhi = pasObj.phi();
    double pasPt  = pasObj.pt();

    double dRval = deltaR(proEta, proPhi, pasEta, pasPhi);
    double dPtRel = 999.0;
    if( proPt > 0.0 ) dPtRel = fabs( pasPt - proPt )/proPt;
    // If we are comparing two objects for which the candidates should
    // be exactly the same, cut hard. Otherwise take cuts from user.
    if( exact ) return ( dRval < 1e-3 && dPtRel < 1e-3 );
    else        return ( dRval < 0.025 && dPtRel < 0.025 );
}

DEFINE_FWK_MODULE(IsolationAnalyzer);
