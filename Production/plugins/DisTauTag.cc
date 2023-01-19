/*
 * Example plugin to demonstrate the direct single-threaded inference with TensorFlow 2.
 */

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
//#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

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

void test_vector(std::vector<float>& values) {
    for (auto& value : values) {
        if (std::isnan(value)) {
            throw runtime_error("DisTauTag score output: NaN detected.");
        } else if (std::isinf(value)) {
            throw runtime_error("DisTauTag score output: Infinity detected.");
        } else if (!std::isfinite(value)) {
            throw runtime_error("DisTauTag score output: Non-standard value detected.");
        }
    }
}

class DisTauTag : public edm::stream::EDProducer<> {
public:
    explicit DisTauTag(const edm::ParameterSet&);
    ~DisTauTag(){};

    // static void fillDescriptions(edm::ConfigurationDescriptions&);

    template<typename Scalar>
    static Scalar getDeltaPhi(Scalar phi1, Scalar phi2);

    static void fill_zero(tensorflow::Tensor&);

private:
    void beginStream(edm::StreamID) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override;

    template <typename FeatureT>
    const float Scale(const Int_t, const Float_t, const bool);
    void saveInputs(const tensorflow::Tensor& tensor, const std::string& block_name);

    std::string graphPath_;

    edm::EDGetTokenT<pat::JetCollection> jets_token;
    edm::EDGetTokenT<pat::PackedCandidateCollection> cands_token;

    const bool save_inputs_;

    tensorflow::GraphDef* graphDef_;
    tensorflow::Session* session_;

    std::ofstream* json_file_;
    
};

// void DisTauTag::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//     // defining this function will lead to a *_cfi file being generated when compiling
//     edm::ParameterSetDescription desc;
//     desc.add<std::string>("graphPath");
//     desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
//     desc.add<edm::InputTag>("pfCandidates", edm::InputTag("packedPFCandidates"));
//     descriptions.addWithDefaultLabel(desc);
// }

DisTauTag::DisTauTag(const edm::ParameterSet& config)
    : graphPath_(config.getParameter<std::string>("graphPath")),
      jets_token(consumes<pat::JetCollection>(config.getParameter<edm::InputTag>("jets"))),
      cands_token(consumes<pat::PackedCandidateCollection>(config.getParameter<edm::InputTag>("pfCandidates"))),
      save_inputs_(config.getParameter<bool>("save_inputs")),
      graphDef_(nullptr),
      session_(nullptr) {
  
  produces<edm::ValueMap<float>>("score0");
  produces<edm::ValueMap<float>>("score1");
  
  tensorflow::setLogging("2");
}

void DisTauTag::beginStream(edm::StreamID) {

    using boost::filesystem::is_regular_file;

    if(!is_regular_file(graphPath_))
        throw std::runtime_error("Error: DisTauTag model can not be opened, path=" + graphPath_);

    // load the graph
    graphDef_ = tensorflow::loadGraphDef(graphPath_);
    // create a new session and add the graphDef
    session_ = tensorflow::createSession(graphDef_);

}

void DisTauTag::endStream() {
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

void DisTauTag::saveInputs(const tensorflow::Tensor& tensor, const std::string& block_name)
{
    int tau_n = tensor.shape().dim_size(0);
    int pf_n = tensor.shape().dim_size(1);
    int ftr_n = tensor.shape().dim_size(2);

    (*json_file_) << "\"" << block_name <<  "\":[";
    for(int tau_idx=0; tau_idx<tau_n; tau_idx++)
    {
        (*json_file_) << "[";
        for(int pf_idx=0; pf_idx<pf_n; pf_idx++)
        {
            (*json_file_) << "[";
            for(int ftr_idx=0; ftr_idx<ftr_n; ftr_idx++)
            {
                (*json_file_) << std::setprecision(7) << std::fixed << tensor.tensor<float, 3>()(tau_idx, pf_idx, ftr_idx);
                if(ftr_idx<ftr_n-1) (*json_file_) << ", ";
            }
            (*json_file_) << "]";
            if(pf_idx<pf_n-1) (*json_file_) << ", ";
        }
        (*json_file_) << "]";
        if(tau_idx<tau_n-1) (*json_file_) << ", ";
    }
    (*json_file_) << "]";
}

void DisTauTag::produce(edm::Event& event, const edm::EventSetup& setup) {

    auto jets = getHandle(event, jets_token);
    auto cands = getHandle(event, cands_token);

    const size_t jets_size = jets->size();
    
    std::vector <Float_t> v_score0(jets_size, -9);
    std::vector <Float_t> v_score1(jets_size, -9);

    // step 1: get jets   
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
            getVecRef(input_1, Br::pfCand_pt                   ,static_cast<Float_t>(daughter->polarP4().pt()));
            getVecRef(input_1, Br::pfCand_eta                  ,static_cast<Float_t>(daughter->polarP4().eta()));
            getVecRef(input_1, Br::pfCand_phi                  ,static_cast<Float_t>(daughter->polarP4().phi()));
            getVecRef(input_1, Br::pfCand_mass                 ,static_cast<Float_t>(daughter->polarP4().mass()));
            getVecRef(input_1, Br::pfCand_charge               ,static_cast<Int_t>(daughter->charge()));
            getVecRef(input_1, Br::pfCand_puppiWeight          ,static_cast<Float_t>(daughter->puppiWeight()));
            getVecRef(input_1, Br::pfCand_puppiWeightNoLep     ,static_cast<Float_t>(daughter->puppiWeightNoLep()));
            getVecRef(input_1, Br::pfCand_lostInnerHits        ,static_cast<Int_t>(daughter->lostInnerHits()));
            getVecRef(input_1, Br::pfCand_nPixelHits           ,static_cast<Int_t>(daughter->numberOfPixelHits()));
            getVecRef(input_1, Br::pfCand_nHits                ,static_cast<Int_t>(daughter->numberOfHits()));
            getVecRef(input_1, Br::pfCand_caloFraction         ,static_cast<Float_t>(daughter->caloFraction()));
            getVecRef(input_1, Br::pfCand_hcalFraction         ,static_cast<Float_t>(daughter->hcalFraction()));
            getVecRef(input_1, Br::pfCand_rawCaloFraction      ,static_cast<Float_t>(daughter->rawCaloFraction()));
            getVecRef(input_1, Br::pfCand_rawHcalFraction      ,static_cast<Float_t>(daughter->rawHcalFraction()));
            
            getVecRef(input_1, Br::pfCand_hasTrackDetails      ,static_cast<Int_t>(daughter->hasTrackDetails()));

            if( daughter->hasTrackDetails() )
            {   
                if(std::isfinite(daughter->dz()))        getVecRef(input_1, Br::pfCand_dz,       static_cast<Float_t>(daughter->dz()));
                if(std::isfinite(daughter->dzError()))   getVecRef(input_1, Br::pfCand_dz_error, static_cast<Float_t>(daughter->dzError()));
                if(std::isfinite(daughter->dxyError()))  getVecRef(input_1, Br::pfCand_dxy_error,static_cast<Float_t>(daughter->dxyError()));

                getVecRef(input_1, Br::pfCand_dxy,        static_cast<Float_t>(daughter->dxy()));
                getVecRef(input_1, Br::pfCand_track_chi2, static_cast<Float_t>(daughter->bestTrack()->chi2()));
                getVecRef(input_1, Br::pfCand_track_ndof, static_cast<Float_t>(daughter->bestTrack()->ndof()));
            }
            
            Float_t jet_eta = jet_p4.eta();
            Float_t jet_phi = jet_p4.phi();
            getVecRef(input_1, PfCand_Features::pfCand_deta, static_cast<Float_t>(daughter->polarP4().eta()) - jet_eta);
            getVecRef(input_1, PfCand_Features::pfCand_dphi, getDeltaPhi<Float_t>(static_cast<Float_t>(daughter->polarP4().phi()), jet_phi));
    
        }

        {   // Categorical features
            typedef PfCandCategorical_Features Br;
            getVecRef(input_2, Br::pfCand_particleType         ,static_cast<Int_t>(TranslatePdgIdToPFParticleType(daughter->pdgId())));
            getVecRef(input_2, Br::pfCand_pvAssociationQuality ,static_cast<Int_t>(daughter->pvAssociationQuality()));
            getVecRef(input_2, Br::pfCand_fromPV               ,static_cast<Int_t>(daughter->fromPV()));
        }

        ++tensor_idx; 

      }

      // define the output and run
      std::vector<tensorflow::Tensor> outputs;
      tensorflow::run(session_,
                    {{"input_1", input_1}, {"input_2", input_2}},
                    {"final_out"}, &outputs);

      // print the output
    //   std::cout << " jet -> " << jetIndex
    //             << " jet_pt -> " << jet_p4.pt()
    //             << " n_pfCand ->" << nDaughters
    //             << " score -> " << outputs[0].flat<float>()(0)
    //             << " " << outputs[0].flat<float>()(1) << std::endl;
    
      v_score0.at(jetIndex) = outputs[0].flat<float>()(0);
      v_score1.at(jetIndex) = outputs[0].flat<float>()(1);

      if (save_inputs_) {

        std::string json_file_name = "distag_"
                + std::to_string(event.id().run()) + "_"
                + std::to_string(event.id().luminosityBlock()) + "_"
                + std::to_string(event.id().event()) + "_" +
                + "jet_" + std::to_string(jetIndex) + ".json";

        json_file_ = new std::ofstream(json_file_name.data());

        (*json_file_) << "{";

        saveInputs(input_1, "PfCand");
        (*json_file_) << ", ";
        saveInputs(input_2, "PfCandCategorical");
        (*json_file_) << ", \"Output\":["
                        << outputs[0].flat<float>()(0) << ","
                        << outputs[0].flat<float>()(1)
                        << "]";

        (*json_file_) << "}";

        delete json_file_;
      }
    }
    
    test_vector(v_score0);
    test_vector(v_score1);
    
    std::unique_ptr<edm::ValueMap<float>> vm_score0(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_score0(*vm_score0);
    filler_score0.insert(jets, v_score0.begin(), v_score0.end());
    filler_score0.fill();
    event.put(std::move(vm_score0), "score0"); // jet probability
    
    
    std::unique_ptr<edm::ValueMap<float>> vm_score1(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler_score1(*vm_score1);
    filler_score1.insert(jets, v_score1.begin(), v_score1.end());
    filler_score1.fill();
    event.put(std::move(vm_score1), "score1"); // tau probability


}

DEFINE_FWK_MODULE(DisTauTag);
