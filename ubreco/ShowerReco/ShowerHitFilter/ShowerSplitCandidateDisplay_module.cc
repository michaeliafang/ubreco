#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include <iostream>
#include <memory>
#include <vector>

// Display-only copy: distinct module label and default-instance associations
// for viewers that do not retain the source product's instance name.
class ShowerSplitCandidateDisplay : public art::EDProducer {
public:
  explicit ShowerSplitCandidateDisplay(const fhicl::ParameterSet& p)
    : EDProducer(p), input_(p.get<art::InputTag>("ClusterProducer")) {
    produces<std::vector<recob::Cluster>>();
    produces<art::Assns<recob::Cluster,recob::Hit>>();
  }
  void produce(art::Event& e) override {
    const auto input = e.getValidHandle<std::vector<recob::Cluster>>(input_);
    const art::FindManyP<recob::Hit> hits(input,e,input_);
    auto clusters = std::make_unique<std::vector<recob::Cluster>>(*input);
    auto links = std::make_unique<art::Assns<recob::Cluster,recob::Hit>>();
    art::PtrMaker<recob::Cluster> makeCluster(e);
    std::size_t count = 0;
    for (std::size_t i=0; i<input->size(); ++i)
      for (const auto& hit : hits.at(i)) {
        links->addSingle(makeCluster(i),hit);
        ++count;
      }
    std::cout << "[ShowerSplitCandidateDisplay] event " << e.id()
      << " clusters=" << clusters->size() << " hit_associations=" << count << std::endl;
    e.put(std::move(clusters));
    e.put(std::move(links));
  }
private:
  art::InputTag input_;
};
DEFINE_ART_MODULE(ShowerSplitCandidateDisplay)
