#include "LLStauProduction/Modules/interface/Utils.h"


namespace Utils
{
    bool isTauSignalCand(const pat::Tau& tau, const pat::PackedCandidate& cand)
    {
        for(const auto& signalCandBase : tau.signalCands())
        {
            auto signalCand = dynamic_cast<const pat::PackedCandidate*>(signalCandBase.get());
            
            
            // For some bloody confounding reason the pointers do not match.
            // So do the computationally more expensive job of checking the dPt and dR.
            // Now it works: READ cand as `const auto &cand = cand_handle->at(icand)` in the code which calls this method.
            // I.e., read as a const reference.
            
            //if(
            //    std::fabs(cand.pt() - signalCandBase->pt()) < 1e-5 &&
            //    reco::deltaR(signalCandBase->eta(), signalCandBase->phi(), cand.eta(), cand.phi()) < 1e-5
            //)
            //{
            //    double dR = reco::deltaR(signalCandBase->eta(), signalCandBase->phi(), cand.eta(), cand.phi());
            //    double dPt = std::fabs(cand.pt() - signalCandBase->pt());
            //    printf(
            //        "Found match: dR(pfcand, tausig) %0.4e, abs(pfcand_pt-tausig_pt) (%4e-%4e) %0.4e, %d \n",
            //        dR, cand.pt(), signalCandBase->pt(), dPt, (int) cand.numberOfSourceCandidatePtrs()
            //    );
            //    
            //    std::cout << signalCand << " " << (&cand) << " " << signalCandBase.get() << " " << (const reco::Candidate*)(&cand) << " " << (const pat::PackedCandidate*)signalCandBase.get() << "\n";
            //    
            //    return true;
            //}
            
            
            if(signalCand == &cand)
            {
                //printf("Pointers match. \n");
                return true;
            }
        }
        
        return false;
    }
    
    bool isTauIsoCand(const pat::Tau& tau, const pat::PackedCandidate& cand)
    {
        for(const auto& isoCandBase : tau.isolationCands())
        {
            //if(
            //    std::fabs(cand.pt() - isoCandBase->pt()) < 1e-5 &&
            //    reco::deltaR(isoCandBase->eta(), isoCandBase->phi(), cand.eta(), cand.phi()) < 1e-5
            //)
            //{
            //    return true;
            //}
            
            auto isoCand = dynamic_cast<const pat::PackedCandidate*>(isoCandBase.get());
            
            if(isoCand == &cand)
            {
                return true;
            }
        }
        
        return false;
    }
    
    bool isTauLeadChHadCand(const pat::Tau& tau, const pat::PackedCandidate& cand)
    {
        //return (
        //    std::fabs(cand.pt() - tau.leadChargedHadrCand()->pt()) < 1e-5 &&
        //    reco::deltaR(tau.leadChargedHadrCand()->eta(), tau.leadChargedHadrCand()->phi(), cand.eta(), cand.phi()) < 1e-5
        //);
        
        auto leadChargedHadrCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
        return leadChargedHadrCand == &cand;
    }
}
