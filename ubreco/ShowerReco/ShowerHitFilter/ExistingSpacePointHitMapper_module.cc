#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ubreco/ShowerReco/ShowerReco3D/Recovery/ExactHitMapping.h"
#include <memory>
#include <chrono>
#include <set>
#include <stdexcept>

namespace {
showerreco::recovery::HitFingerprint fingerprint(const recob::Hit& h) {
  const auto& w = h.WireID();
  return {{double(h.Channel()),double(w.Cryostat),double(w.TPC),double(w.Plane),double(w.Wire),
    double(h.View()),double(h.SignalType()),double(h.StartTick()),double(h.EndTick()),
    h.PeakTime(),h.SigmaPeakTime(),h.RMS(),h.PeakAmplitude(),h.SigmaPeakAmplitude(),
    h.Integral(),h.SigmaIntegral(),h.ROISummedADC(),double(h.Multiplicity()),
    double(h.LocalIndex()),h.GoodnessOfFit(),double(h.DegreesOfFreedom())}};
}
}

// Reuse persisted spacepoints without copying or refitting them. Exact matching
// is appropriate because both intervening selectors copy recob::Hit unchanged.
class ExistingSpacePointHitMapper : public art::EDProducer {
public:
  explicit ExistingSpacePointHitMapper(const fhicl::ParameterSet& p)
    : EDProducer(p), points_(p.get<art::InputTag>("SpacePointProducer")),
      source_(p.get<art::InputTag>("SourceHitProducer")),
      selected_(p.get<art::InputTag>("SelectedHitProducer")) {
    produces<art::Assns<recob::SpacePoint,recob::Hit>>();
    produces<std::vector<art::Ptr<recob::Hit>>>("source");
  }
  void produce(art::Event& e) override {
    const auto started = std::chrono::steady_clock::now();
    using namespace showerreco::recovery;
    const auto points = e.getValidHandle<std::vector<recob::SpacePoint>>(points_);
    const auto source = e.getValidHandle<std::vector<recob::Hit>>(source_);
    const auto selected = e.getValidHandle<std::vector<recob::Hit>>(selected_);
    const art::FindManyP<recob::Hit> links(points,e,points_);
    std::vector<HitFingerprint> src, dst;
    for (const auto& h : *source) src.push_back(fingerprint(h));
    for (const auto& h : *selected) dst.push_back(fingerprint(h));
    const auto mapping = mapCopiedHits(src,dst);
    auto output = std::make_unique<art::Assns<recob::SpacePoint,recob::Hit>>();
    // Parallel to selected hits; null entries mark unresolved correspondence.
    auto provenance = std::make_unique<std::vector<art::Ptr<recob::Hit>>>(selected->size());
    for (std::size_t i=0; i<mapping.sourceToSelected.size(); ++i)
      if (mapping.sourceToSelected[i]!=noHit)
        provenance->at(mapping.sourceToSelected[i]) = art::Ptr<recob::Hit>(source,i);
    std::size_t kept = 0, associations = 0;
    for (std::size_t i=0; i<points->size(); ++i) {
      std::set<std::size_t> target;
      for (const auto& h : links.at(i)) {
        if (h.id()!=source.id())
          throw std::runtime_error("ExistingSpacePointHitMapper: sps hit ProductID differs from SourceHitProducer");
        const auto j = mapping.sourceToSelected.at(h.key());
        if (j!=noHit) target.insert(j);
      }
      kept += !target.empty();
      for (auto j : target) {
        output->addSingle(art::Ptr<recob::SpacePoint>(points,i),art::Ptr<recob::Hit>(selected,j));
        ++associations;
      }
    }
    mf::LogInfo("ExistingSpacePointHitMapper") << "event " << e.id()
      << " source_hits=" << source->size() << " selected_hits=" << selected->size()
      << " matched=" << mapping.matched << " missing_selected=" << mapping.missing
      << " ambiguous_selected=" << mapping.ambiguous << " invalid_selected=" << mapping.invalid
      << " mapped_points=" << kept << " associations=" << associations
      << " mapping_ms=" << std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
    e.put(std::move(output));
    e.put(std::move(provenance),"source");
  }
private:
  art::InputTag points_, source_, selected_;
};
DEFINE_ART_MODULE(ExistingSpacePointHitMapper)
