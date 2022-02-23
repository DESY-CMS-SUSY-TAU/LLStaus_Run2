#ifndef Utils_H
#define Utils_H


// system include files
# include <algorithm>
# include <memory>

// user include files
# include "CommonTools/UtilAlgos/interface/TFileService.h"
# include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
# include "DataFormats/Common/interface/MapOfVectors.h"
# include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
# include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
# include "DataFormats/EgammaCandidates/interface/Photon.h"
# include "DataFormats/FWLite/interface/ESHandle.h"
# include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
# include "DataFormats/HepMCCandidate/interface/GenParticle.h"
# include "DataFormats/JetReco/interface/PFJet.h"
# include "DataFormats/Math/interface/LorentzVector.h"
# include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
# include "DataFormats/PatCandidates/interface/Jet.h"
# include "DataFormats/PatCandidates/interface/MET.h"
# include "DataFormats/PatCandidates/interface/Muon.h"
# include "DataFormats/PatCandidates/interface/Electron.h"
# include "DataFormats/PatCandidates/interface/Photon.h"
# include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
# include "DataFormats/PatCandidates/interface/Tau.h"
# include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
# include "DataFormats/TrackReco/interface/Track.h"
# include "DataFormats/TrackReco/interface/TrackFwd.h"
# include "DataFormats/VertexReco/interface/Vertex.h"
# include "FWCore/Framework/interface/ConsumesCollector.h"
# include "FWCore/Framework/interface/Event.h"
# include "FWCore/Framework/interface/ESHandle.h"
# include "FWCore/Framework/interface/Frameworkfwd.h"
# include "FWCore/Framework/interface/MakerMacros.h"
# include "FWCore/Framework/interface/one/EDAnalyzer.h"
# include "FWCore/ParameterSet/interface/ParameterSet.h"
# include "FWCore/ServiceRegistry/interface/Service.h"
# include "FWCore/Utilities/interface/InputTag.h"
# include "Geometry/Records/interface/CaloGeometryRecord.h"
# include "Geometry/Records/interface/IdealGeometryRecord.h"
# include "RecoEgamma/EgammaTools/interface/MVAVariableHelper.h"
# include "RecoEgamma/EgammaTools/interface/MVAVariableManager.h"
# include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
# include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
# include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
# include "SimDataFormats/CaloHit/interface/PCaloHit.h"
# include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
# include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

# include <Compression.h>
# include <Math/VectorUtil.h>
# include <TH1F.h>
# include <TH2F.h>
# include <TLorentzVector.h>
# include <TMatrixD.h>
# include <TTree.h>
# include <TVector2.h>
# include <TVectorD.h>


# include <CLHEP/Matrix/Matrix.h>
# include <CLHEP/Vector/LorentzVector.h>
# include <CLHEP/Vector/ThreeVector.h>
# include <CLHEP/Vector/ThreeVector.h>


namespace Utils
{
    bool isTauSignalCand(const pat::Tau& tau, const pat::PackedCandidate& cand);
    bool isTauIsoCand(const pat::Tau& tau, const pat::PackedCandidate& cand);
    bool isTauLeadChHadCand(const pat::Tau& tau, const pat::PackedCandidate& cand);
}


#endif
