// -*- C++ -*-
//
// Package:    Exclusivity/ExlAnalyzer
// Class:      ExlAnalyzer
// 
/**\class ExlAnalyzer ExlAnalyzer.cc Exclusivity/ExlAnalyzer/plugins/ExlAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Georgios Krintiras
//         Created:  Fri, 21 Aug 2020 20:29:52 GMT
//
//


// system include files
#include <memory>
#include <ostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ExlAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ExlAnalyzer(const edm::ParameterSet&);
      ~ExlAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<edm::SortedCollection<CaloTower>>        CaloTowerCollection_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsCollection_;
  edm::EDGetTokenT<edm::View<reco::Photon>> photonsCollection_;
  //enum EDet { kEB, kEE, kHB, kHE, kHFm, kHFp, nDets };
  std::vector<std::string> detNames = { "EB", "EE", "HB", "HE", "HFm", "HFp", "unknown" };

  std::vector<std::pair<double,double>> detLimits = {
    {0    , 1.4442 }, // EB (Exclude transition region between calo barrel and endcap)
    {1.566, 3.0   }, // EE (Exclude transition region between calo barrel and endcap)
    {0    , 1.3   }, // HB
    {1.3  , 3.0   }, // HE
    {-5.2 ,-3.0   }, // HFm
    {3.0  , 5.2   }, // HFp*/
  };

  enum ECaloType { kEB, kEE, kHB, kHE, kHFp, kHFm, nCaloTypes };
  //constexpr initializer_list<ECaloType> calotypes = {kEB, kEE, kHB, kHE, kHFp, kHFm};
  ECaloType GetTowerSubdetHad(double&) const;
  ECaloType GetTowerSubdetEm(double&)  const;

  const std::map<ECaloType, std::string> caloName = {
    { kEB  , "EB"  },
    { kEE  , "EE"  },
    { kHB  , "HB"  },
    { kHE  , "HE"  },
    { kHFp , "HFp" },
    { kHFm , "HFm" },
  };

  const std::map<ECaloType, double> noiseThreshold = {
    { kEB  , 0.7 },
    { kEE  , 7.5 },
    { kHB  , 2.8  },
    { kHE  , 2.4  },
    { kHFp , 7.2  },
    { kHFm , 7.5  },
  };

  const double maxEtaEB = detLimits.at(0).second;
  const double minEtaEE = detLimits.at(1).first;
  const double maxEtaEE = detLimits.at(1).second;

  const double maxEtaHB = detLimits.at(2).second;
  const double minEtaHE = detLimits.at(3).first;
  const double maxEtaHE = detLimits.at(3).second;

  const double minEtaHF = abs(detLimits.at(4).first);
  const double maxEtaHF = abs(detLimits.at(4).second);

  bool IsInHEM(double&eta, double&phi);
  bool IsOverlappingWithElectron(const edm::Event& , ECaloType , double&, double&);
  bool IsOverlappingWithPhoton(const edm::Event& , ECaloType , double&, double&);

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ExlAnalyzer::ExlAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  //CaloTowerCollection_     = consumes<edm::SortedCollection<CaloTower>>(iConfig.getParameter<edm::InputTag>("recoCaloTower"));
  CaloTowerCollection_     = consumes<edm::SortedCollection<CaloTower>>(edm::InputTag("towerMaker"));
  electronsCollection_ = consumes<edm::View<reco::GsfElectron>>(edm::InputTag("slimmedElectrons"));
  photonsCollection_ = consumes<edm::View<reco::Photon>>(edm::InputTag("slimmedPhotons"));
  usesResource("TFileService");
   
}


ExlAnalyzer::~ExlAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExlAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
  //Calo towers
   edm::Handle<edm::SortedCollection<CaloTower>> CaloTowerHandle;
   iEvent.getByToken(CaloTowerCollection_, CaloTowerHandle);

   double eta;
   double phi;
   double energy;
   double energyEm;
   double energyHad;
   //double energyTransverse;
   
   for (edm::SortedCollection<CaloTower>::const_iterator calo = CaloTowerHandle->begin(); calo != CaloTowerHandle->end(); ++calo) {
     
     energy=calo->energy();
     energyEm=calo->emEnergy();
     energyHad= calo->hadEnergy();
     //energyTransverse=calo->et();
     phi=calo->phi();
     eta=calo->eta();

     ECaloType subdetHad = GetTowerSubdetHad(eta);
     if(subdetHad==kHFp || subdetHad==kHFm){ // Check HCAL and HF exclusivity
       if(energy > noiseThreshold.at(subdetHad)){ 
	 std::cout<< " rejected due to "<<caloName.at(subdetHad)<<std::endl;
       }
     }
     
     if(subdetHad==kHB || subdetHad==kHE){ // Check HCAL and HF exclusivity                                                                                                                    
       if(energyHad > noiseThreshold.at(subdetHad)){
	 std::cout<< " rejected due to "<<caloName.at(subdetHad)<<std::endl;
       }
     }


     ECaloType subdetEm = GetTowerSubdetEm(eta);
     if(subdetEm==kEB){ // Check EB exclusivity 
       if(IsOverlappingWithElectron(iEvent, subdetEm, eta, phi)) continue; //Prob a good candidate
       if(IsOverlappingWithPhoton(iEvent, subdetEm, eta, phi)) continue; //Prob a good candidate 
       if(energyEm > noiseThreshold.at(subdetEm)){
	 std::cout<< " rejected due to "<<caloName.at(subdetEm)<<std::endl;
       }
     }
     if(subdetEm==kEE){ // Check EE exclusivity
       if(fabs(eta) >  2.3) continue; //Don't look at towers that are in very noisy region of EE
       if(IsInHEM(eta,phi)) continue; // Don't look at towers that are in noisy HEM region (only for 2018!!)
       if(IsOverlappingWithElectron(iEvent, subdetEm, eta, phi)) continue; //Prob a good candidate                                                                                                     
       if(IsOverlappingWithPhoton(iEvent, subdetEm, eta, phi)) continue; //Prob a good candidate   

       if(energyEm > noiseThreshold.at(subdetEm)){
	 std::cout<< " rejected due to "<<caloName.at(subdetEm)<<std::endl;
       }
     }
   }// calo tower loop


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

ExlAnalyzer::ECaloType 
ExlAnalyzer::GetTowerSubdetHad(double&eta) const
{
  if(eta > -maxEtaHF && eta < -minEtaHF) return kHFm;
  if(eta >  minEtaHF && eta <  maxEtaHF) return kHFp;
  if(fabs(eta) > 0 && fabs(eta) < maxEtaHB) return kHB;
  if(fabs(eta) > minEtaHE && fabs(eta) < maxEtaHE) return kHE;

  return nCaloTypes;
}

ExlAnalyzer::ECaloType 
ExlAnalyzer::GetTowerSubdetEm(double&eta) const
{
  if(eta > -maxEtaHF && eta < -minEtaHF) return kHFm;
  if(eta >  minEtaHF && eta <  maxEtaHF) return kHFp;
  
  if(fabs(eta) > 0 && fabs(eta) < maxEtaEB) return kEB;
  if(fabs(eta) > minEtaEE && fabs(eta) < maxEtaEE) return kEE;
  
  return nCaloTypes;
}

bool
ExlAnalyzer::IsInHEM(double&eta, double&phi)
{
  if(eta < -1.39 &&
     eta > -maxEtaEE &&
     phi > -1.6 &&
     phi < -0.9) return true;
  
  return false;
}

bool 
ExlAnalyzer::IsOverlappingWithElectron(const edm::Event& iEvent, ECaloType subdet, double&eta, double&phi)
{

  //Electrons
  edm::Handle<edm::View<reco::GsfElectron>> gsfElectrons;
  iEvent.getByToken(electronsCollection_, gsfElectrons);

  // loop over electrons
  bool overlapsWithElectron = false;
  double maxDeltaEta = (subdet == kEB ) ? 0.15 : 0.15;
  double maxDeltaPhi = (subdet == kEB ) ? 0.7 : 0.4;

  for (auto ele = gsfElectrons->begin(); ele != gsfElectrons->end(); ++ele) {

    if(ele->pt() <  2.0) continue;

    double deltaEta = fabs(ele->superCluster()->eta() - eta);
    double deltaPhi = fabs(ele->superCluster()->phi() - phi);
    
    if(deltaEta < maxDeltaEta && deltaPhi < maxDeltaPhi){
      overlapsWithElectron = true;
      break;
    }
  }
  return overlapsWithElectron;
}

bool
ExlAnalyzer::IsOverlappingWithPhoton(const edm::Event& iEvent, ECaloType subdet, double&eta, double&phi)
{

  //Photons
  edm::Handle<edm::View<reco::Photon>> recoPhotons;
  iEvent.getByToken(photonsCollection_, recoPhotons);

  // loop over photons                                                                                                                                                                                    
  bool overlapsWithPhoton = false;
  double maxDeltaEta = (subdet == kEB ) ? 0.15 : 0.15;
  double maxDeltaPhi = (subdet == kEB ) ? 0.7 : 0.4;

  for (auto pho = recoPhotons->begin(); pho != recoPhotons->end(); ++pho) {

    if(pho->et() <  2.0) continue;

    double deltaEta = fabs(pho->superCluster()->eta() - eta);
    double deltaPhi = fabs(pho->superCluster()->phi() - phi);

    if(deltaEta < maxDeltaEta && deltaPhi < maxDeltaPhi){
      overlapsWithPhoton = true;
      break;
    }
  }
  return overlapsWithPhoton;
}

// ------------ method called once each job just before starting event loop  ------------
void 
ExlAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExlAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExlAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExlAnalyzer);
