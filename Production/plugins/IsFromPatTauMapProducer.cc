// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      IsFromPatTauMapProducer
//
/**\class IsFromPatTauMapProducer IsFromPatTauMapProducer.cc PhysicsTools/NanoAOD/plugins/IsFromPatTauMapProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maria Giulia Ratti (ETHZ) [mratti]
//         Created:  Thu, 22 Nov 2018 12:34:48 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "LLStaus_Run2/Production/interface/Utils.h"

//
// class declaration
//

class IsFromPatTauMapProducer : public edm::global::EDProducer<>
{
    public:
        explicit IsFromPatTauMapProducer(const edm::ParameterSet& iConfig):
            pc_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCandidates"))), // pf candidates
            tt_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("patTaus"))) // taus
        {
            produces<edm::ValueMap<int>>("isTauIdxSignalCand"); // name of the value map that I want to actually produce
            produces<edm::ValueMap<int>>("isTauIdxIsoCand"); // name of the value map that I want to actually produce
            produces<edm::ValueMap<int>>("isTauIdxLeadChHadCand"); // name of the value map that I want to actually produce
        }
        ~IsFromPatTauMapProducer() override {};
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        
    private:
        void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
        
        
        // ----------member data ---------------------------
        edm::EDGetTokenT<pat::PackedCandidateCollection> pc_;
        edm::EDGetTokenT<pat::TauCollection> tt_;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// member functions
//

// ------------ method called to produce the data  ------------
void IsFromPatTauMapProducer::produce(edm::StreamID streamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    // packedPFCandidate collection
    edm::Handle<pat::PackedCandidateCollection> pc_handle;
    iEvent.getByToken( pc_, pc_handle );
    
    // tau collection
    edm::Handle<pat::TauCollection> tau_handle;
    iEvent.getByToken( tt_, tau_handle );
    
    // the map cannot be filled straight away, so create an intermediate vector
    unsigned int Npc = pc_handle->size();
    std::vector<int> v_isTauIdxSignalCand(Npc, -1);
    std::vector<int> v_isTauIdxIsoCand(Npc, -1);
    std::vector<int> v_isTauIdxLeadChHadCand(Npc, -1);
    
    unsigned int Ntau = tau_handle->size();
    
    for (unsigned int ipc=0; ipc<Npc; ipc++)
    {
        const auto &pc = pc_handle->at(ipc);
        
        int isTauIdxSignalCand = -1;
        int isTauIdxIsoCand = -1;
        int isTauIdxLeadChHadCand = -1;
        
        for (unsigned int itau=0; itau<Ntau; itau++)
        {
            const auto &tau = tau_handle->at(itau);
            
            // Only check if not already matched
            if (isTauIdxSignalCand < 0 && Utils::isTauSignalCand(tau, pc))
            {
                isTauIdxSignalCand = itau;
            }
            
            if (isTauIdxIsoCand < 0 && Utils::isTauIsoCand(tau, pc))
            {
                isTauIdxIsoCand = itau;
            }
            
            if (isTauIdxLeadChHadCand < 0 && Utils::isTauLeadChHadCand(tau, pc))
            {
                isTauIdxLeadChHadCand = itau;
            }
        }
        
        //if (isTauIdxSignalCand >= 0 || isTauIdxIsoCand >= 0 || isTauIdxLeadChHadCand >= 0)
        //{
        //    printf("Ntau %d, ipc %d, isTauIdxSignalCand %d, isTauIdxIsoCand %d, isTauIdxLeadChHadCand %d \n", (int) Ntau, (int) ipc, isTauIdxSignalCand, isTauIdxIsoCand, isTauIdxLeadChHadCand);
        //}
        
        v_isTauIdxSignalCand[ipc] = isTauIdxSignalCand;
        v_isTauIdxIsoCand[ipc] = isTauIdxIsoCand;
        v_isTauIdxLeadChHadCand[ipc] = isTauIdxLeadChHadCand;
    }
    
    
    std::unique_ptr<edm::ValueMap<int>> vm_isTauIdxSignalCand(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler filler_isTauIdxSignalCand(*vm_isTauIdxSignalCand);
    filler_isTauIdxSignalCand.insert(pc_handle, v_isTauIdxSignalCand.begin(), v_isTauIdxSignalCand.end());
    filler_isTauIdxSignalCand.fill();
    iEvent.put(std::move(vm_isTauIdxSignalCand), "isTauIdxSignalCand");
    
    
    std::unique_ptr<edm::ValueMap<int>> vm_isTauIdxIsoCand(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler filler_isTauIdxIsoCand(*vm_isTauIdxIsoCand);
    filler_isTauIdxIsoCand.insert(pc_handle, v_isTauIdxIsoCand.begin(), v_isTauIdxIsoCand.end());
    filler_isTauIdxIsoCand.fill();
    iEvent.put(std::move(vm_isTauIdxIsoCand), "isTauIdxIsoCand");
    
    
    std::unique_ptr<edm::ValueMap<int>> vm_isTauIdxLeadChHadCand(new edm::ValueMap<int>());
    edm::ValueMap<int>::Filler filler_isTauIdxLeadChHadCand(*vm_isTauIdxLeadChHadCand);
    filler_isTauIdxLeadChHadCand.insert(pc_handle, v_isTauIdxLeadChHadCand.begin(), v_isTauIdxLeadChHadCand.end());
    filler_isTauIdxLeadChHadCand.fill();
    iEvent.put(std::move(vm_isTauIdxLeadChHadCand), "isTauIdxLeadChHadCand");
    
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void IsFromPatTauMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{

    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("packedPFCandidates")->setComment("packed PF Candidates collection ");
    desc.add<edm::InputTag>("patTaus")->setComment("tau collection");
    
    descriptions.addWithDefaultLabel(desc);
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsFromPatTauMapProducer);
