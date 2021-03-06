/*
 * Example plugin to demonstrate the direct single-threaded inference with TensorFlow 2.
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "LLStaus_Run2/Production/interface/DisTauTag_input.h"

namespace {
    template <typename T>
    edm::Handle<T> getHandle(const edm::Event& event, const edm::EDGetTokenT<T>& token, bool get = true) {
        edm::Handle<T> handle;
        if(get)
            event.getByToken(token, handle);
        return handle;
    }
}

class DisTauTag : public edm::one::EDAnalyzer<> {
public:
    explicit DisTauTag(const edm::ParameterSet&);
    ~DisTauTag(){};

    static void fillDescriptions(edm::ConfigurationDescriptions&);

    template<typename Scalar>
    static Scalar getDeltaPhi(Scalar phi1, Scalar phi2);

    static void fill_zero(tensorflow::Tensor&);

private:
    void beginJob();
    void analyze(const edm::Event&, const edm::EventSetup&);
    void endJob();

    template <typename FeatureT>
    const float Scale(const Int_t, const Float_t, const bool);

    std::string graphPath_;

    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::PackedCandidateCollection> cands_token;

    tensorflow::GraphDef* graphDef_;
    tensorflow::Session* session_;

};

void DisTauTag::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    // defining this function will lead to a *_cfi file being generated when compiling
    edm::ParameterSetDescription desc;
    desc.add<std::string>("graphPath");
    desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
    desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
    descriptions.addWithDefaultLabel(desc);
}

DisTauTag::DisTauTag(const edm::ParameterSet& config)
    : graphPath_(config.getParameter<std::string>("graphPath")),
      jets_token(consumes<pat::JetCollection>(config.getParameter<edm::InputTag>("jets"))),
      cands_token(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("pfCandidates"))),
      graphDef_(nullptr),
      session_(nullptr) {

  tensorflow::setLogging("2");
}

void DisTauTag::beginJob() {

    using boost::filesystem::is_regular_file;

    if(!is_regular_file(graphPath_))
        throw std::runtime_error("Error: DisTauTag model can not be opened, path=" + graphPath_);

    // load the graph
    graphDef_ = tensorflow::loadGraphDef(graphPath_);
    // create a new session and add the graphDef
    session_ = tensorflow::createSession(graphDef_);

}

void DisTauTag::endJob() {
    // close the session
    tensorflow::closeSession(session_);

    // delete the graph
    delete graphDef_;
    graphDef_ = nullptr;
}

template<typename Scalar>
Scalar DisTauTag::getDeltaPhi(Scalar phi1, Scalar phi2)
{
    static constexpr Scalar pi = boost::math::constants::pi<Scalar>();
    Scalar dphi = phi1 - phi2;
    if(dphi > pi)
        dphi -= 2*pi;
    else if(dphi <= -pi)
        dphi += 2*pi;
    return dphi;
};

void DisTauTag::fill_zero(tensorflow::Tensor& tensor)
{
    size_t size_ = 1;
    int num_dimensions = tensor.shape().dims();
    for(int ii_dim=0; ii_dim<num_dimensions; ii_dim++)
        size_=size_*tensor.shape().dim_size(ii_dim);
    
    for(size_t ii=0; ii<size_; ii++)
        tensor.flat<float>()(ii) = 0.0;
}

template <typename FeatureT>
const float DisTauTag::Scale(const Int_t idx, const Float_t value, const bool inner)
{
    return std::clamp((value - FeatureT::mean.at(idx).at(inner)) / FeatureT::std.at(idx).at(inner),
                        FeatureT::lim_min.at(idx).at(inner), FeatureT::lim_max.at(idx).at(inner));
}

void DisTauTag::analyze(const edm::Event& event, const edm::EventSetup& setup) {

    auto jets = getHandle(event, jets_token);
    auto cands = getHandle(event, cands_token);

    // step 1: get jets
    const size_t jets_size = jets->size();
    for(size_t jetIndex = 0; jetIndex < jets_size; ++jetIndex)
    {
      const auto& jet = jets->at(jetIndex);
      const auto& jet_p4 = jet.polarP4();

       // step 2: get jet dughters
      const size_t nDaughters = jet.numberOfDaughters();

      // step 3: sort by pt
      std::vector<size_t> indices(nDaughters);
      std::iota(indices.begin(), indices.end(), 0);
      std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        const auto& daughter_1 = jet.daughterPtr(a);
        const auto& daughter_2 = jet.daughterPtr(b);
        return daughter_1->polarP4().pt() > daughter_2->polarP4().pt();
      });

      // step 4: mapping function for the scaling
      tensorflow::Tensor input_1(tensorflow::DT_FLOAT,
                                tensorflow::TensorShape({1, Setup::nSeq_PfCand, Setup::n_PfCand}));
      tensorflow::Tensor input_2(tensorflow::DT_FLOAT,
                                tensorflow::TensorShape({1, Setup::nSeq_PfCand, Setup::n_PfCandCategorical}));
      fill_zero(input_1);
      fill_zero(input_2);

      size_t daughter_idx = 0;
      size_t tensor_idx = 0;

      while( tensor_idx < Setup::nSeq_PfCand && daughter_idx < nDaughters)
      {
        const auto& jet_daughter = jet.daughterPtr(indices.at(daughter_idx));
        const auto daughter = dynamic_cast<const pat::PackedCandidate*>(jet_daughter.get());
        ++daughter_idx;

        auto getVecRef = [&](tensorflow::Tensor& tensor, auto _fe, Float_t value){
          const int _feature_idx = static_cast<int>(_fe);
          if(_feature_idx < 0) return;
          tensor.tensor<float, 3>()(0, tensor_idx, _feature_idx)
              = Scale<typename  FeaturesHelper<decltype(_fe)>::scaler_type>(_feature_idx, value, false);
        };

        {   // General features
            typedef PfCand_Features Br;
            getVecRef(input_1, Br::pfCand_valid                ,1.0);
            getVecRef(input_1, Br::pfCand_pt                   ,static_cast<float>(daughter->polarP4().pt()));
            getVecRef(input_1, Br::pfCand_eta                  ,static_cast<float>(daughter->polarP4().eta()));
            getVecRef(input_1, Br::pfCand_phi                  ,static_cast<float>(daughter->polarP4().phi()));
            getVecRef(input_1, Br::pfCand_mass                 ,static_cast<float>(daughter->polarP4().mass()));
            getVecRef(input_1, Br::pfCand_charge               ,static_cast<float>(daughter->charge()));
            getVecRef(input_1, Br::pfCand_puppiWeight          ,static_cast<float>(daughter->puppiWeight()));
            getVecRef(input_1, Br::pfCand_puppiWeightNoLep     ,static_cast<float>(daughter->puppiWeightNoLep()));
            getVecRef(input_1, Br::pfCand_lostInnerHits        ,static_cast<float>(daughter->lostInnerHits()));
            getVecRef(input_1, Br::pfCand_nPixelHits           ,static_cast<float>(daughter->numberOfPixelHits()));
            getVecRef(input_1, Br::pfCand_nHits                ,static_cast<float>(daughter->numberOfHits()));
            getVecRef(input_1, Br::pfCand_caloFraction         ,static_cast<float>(daughter->caloFraction()));
            getVecRef(input_1, Br::pfCand_hcalFraction         ,static_cast<float>(daughter->hcalFraction()));
            getVecRef(input_1, Br::pfCand_rawCaloFraction      ,static_cast<float>(daughter->rawCaloFraction()));
            getVecRef(input_1, Br::pfCand_rawHcalFraction      ,static_cast<float>(daughter->rawHcalFraction()));
            
            getVecRef(input_1, Br::pfCand_hasTrackDetails      ,static_cast<float>(daughter->hasTrackDetails()));

            if( daughter->hasTrackDetails() )
            {   
                if(std::isfinite(daughter->dz()))        getVecRef(input_1, Br::pfCand_dz,       static_cast<float>(daughter->dz()));
                if(std::isfinite(daughter->dzError()))   getVecRef(input_1, Br::pfCand_dz_error, static_cast<float>(daughter->dzError()));
                if(std::isfinite(daughter->dxyError()))  getVecRef(input_1, Br::pfCand_dxy_error,static_cast<float>(daughter->dxyError()));

                getVecRef(input_1, Br::pfCand_dxy,        static_cast<float>(daughter->dxy()));
                getVecRef(input_1, Br::pfCand_track_chi2, static_cast<float>(daughter->bestTrack()->chi2()));
                getVecRef(input_1, Br::pfCand_track_ndof, static_cast<float>(daughter->bestTrack()->ndof()));
            }
            
            float jet_eta = jet_p4.eta();
            float jet_phi = jet_p4.phi();
            getVecRef(input_1, PfCand_Features::pfCand_deta, static_cast<float>(daughter->polarP4().eta()) - jet_eta);
            getVecRef(input_1, PfCand_Features::pfCand_dphi, getDeltaPhi<float>(static_cast<float>(daughter->polarP4().eta()), jet_phi));
    
        }

        {   // Categorical features
            typedef PfCandCategorical_Features Br;
            getVecRef(input_2, Br::pfCand_particleType         ,static_cast<float>(TranslatePdgIdToPFParticleType(daughter->pdgId())));
            getVecRef(input_2, Br::pfCand_pvAssociationQuality ,static_cast<float>(daughter->pvAssociationQuality()));
            getVecRef(input_2, Br::pfCand_fromPV               ,static_cast<float>(daughter->fromPV()));
        }

        ++tensor_idx; 

      }

      // define the output and run
      std::vector<tensorflow::Tensor> outputs;
      tensorflow::run(session_,
                    {{"input_1", input_1}, {"input_2", input_2}},
                    {"final_out"}, &outputs);

      // print the output
      std::cout << " jet -> " << jetIndex
                << " n_pfCand ->" << nDaughters
                << " score -> " << outputs[0].flat<float>()(0)
                << " " << outputs[0].flat<float>()(1) << std::endl;
    }

}

DEFINE_FWK_MODULE(DisTauTag);