#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

using namespace std;

class TTLLGenCategoryFilter : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TTLLGenCategoryFilter(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  typedef std::vector<int> vint;
  enum InputType { IN_PartonTop, IN_PseudoTop, IN_Hadron } inputType_;
  enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON } decayMode_;
  const bool doInvert_;

  // For the parton top
  edm::EDGetTokenT<reco::GenParticleCollection> parton_srcToken_;
  edm::EDGetTokenT<int> parton_channelToken_;
  edm::EDGetTokenT<vint> parton_modesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> parton_jetToken_;
  bool vetoTau_;

  // FIXME : classification by pseudotop, hadrons to be added
};

TTLLGenCategoryFilter::TTLLGenCategoryFilter(const edm::ParameterSet& pset):
  doInvert_(pset.getParameter<bool>("invert"))
{
  const auto inputType = pset.getParameter<string>("inputType");
  if ( inputType == "PartonTop" ) inputType_ = IN_PartonTop;
  //else if ( inputType == "PseudoTop" ) inputType_ = IN_PseudoTop;
  //else if ( inputType == "Hadron" ) inputType_ = IN_Hadron;
  else edm::LogError("TTLLGenCategoryFilter") << "Wrong input inputType. Choose among (PartonTop,)";

  const auto addJetType = pset.getParameter<string>("addJetType");

  switch ( inputType_ )
  {
    case IN_PartonTop: {
      vetoTau_ = pset.getParameter<bool>("vetoTau");
      const auto label = pset.getParameter<edm::InputTag>("src");
      const auto labelName = label.label();
      parton_srcToken_ = consumes<reco::GenParticleCollection>(label);
      parton_channelToken_ = consumes<int>(edm::InputTag(labelName, "channel"));
      parton_modesToken_ = consumes<vint>(edm::InputTag(labelName, "modes"));
      parton_jetToken_ = consumes<reco::GenJetCollection>(edm::InputTag(labelName, "qcdJets"));
      break;
    }
    case IN_PseudoTop: {
      break;
    }
    case IN_Hadron: {
      break;
    }
  }

}

bool TTLLGenCategoryFilter::filter(edm::Event& event, const edm::EventSetup&)
{
  bool accept = false;
  switch ( inputType_ )
  {
    case IN_PartonTop: {
      edm::Handle<reco::GenParticleCollection> srcHandle;
      edm::Handle<int> channelHandle;
      edm::Handle<vint> modesHandle;
      edm::Handle<reco::GenJetCollection> jetHandle;

      event.getByToken(parton_srcToken_, srcHandle);
      event.getByToken(parton_channelToken_, channelHandle);
      event.getByToken(parton_modesToken_, modesHandle);
      event.getByToken(parton_jetToken_, jetHandle);

      const int channel = *channelHandle;
      const int mode1 = modesHandle->at(0);
      const int mode2 = modesHandle->at(1);

      // Dilepton channel only
      if ( channel != 3 ) break;

      // Veto tau channel if set
      if ( vetoTau_ )
      {
        if ( mode1 >= CH_TAU_HADRON ) break;
        if ( mode2 >= CH_TAU_HADRON ) break;
      }
    }
    case IN_PseudoTop: {
      break;
    }
    case IN_Hadron: {
      break;
    }
  }

  if ( doInvert_ ) return !accept;
  return accept;
}

