#ifndef UBRECO_SINGLESHOWERRECOVERY_H
#define UBRECO_SINGLESHOWERRECOVERY_H

#include "VertexSplit.h"
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/ClusterMaker.h"
#include "ubreco/ShowerReco/ClusterMerging/CMToolBase/CMatchManager.h"
#include "ubreco/ShowerReco/ShowerReco3D/Base/ShrRecoManager.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include <set>
#include <string>

namespace showerreco {
namespace recovery {

using HitList = std::vector<art::Ptr<recob::Hit>>;
struct RecoveryResult {
  std::vector<recob::Cluster> clusters;
  std::vector<HitList> clusterHits;
  std::vector<std::vector<unsigned int>> matches;
  std::vector<Shower_t> showers;
};

// Mirrors the pre-publication validity checks in ShrReco3D::SaveShower.
inline bool usableShower(const Shower_t& s) {
  return s.fPassedReconstruction && s.fXYZStart.Mag2() != 0. && s.fDCosStart.Mag2() != 0.;
}

class SingleShowerRecovery {
public:
  explicit SingleShowerRecovery(const fhicl::ParameterSet& p) : matcher_(3) {
    cfg_.minHits = p.get<unsigned int>("MinHits", 10);
    cfg_.minUpstreamHits = p.get<unsigned int>("MinUpstreamHits", 3);
    cfg_.maxHits = p.get<unsigned int>("MaxHits", 6000);
    cfg_.neighborRadius = p.get<double>("NeighborRadius", 0.8);
    cfg_.innerRadius = p.get<double>("InnerRadius", 5.);
    cfg_.upstreamRadius = p.get<double>("UpstreamRadius", 20.);
    cfg_.angleStepDegrees = p.get<double>("AngleStepDegrees", 2.);
    cfg_.minOpeningDegrees = p.get<double>("MinOpeningDegrees", 10.);
    cfg_.minChargeFraction = p.get<double>("MinChargeFraction", 0.05);
    cfg_.maxCutFraction = p.get<double>("MaxCutFraction", 0.15);
    cfg_.validate();
    maxCandidates_ = p.get<unsigned int>("MaxCandidatesPerPlane", 2);
    minCoverage_ = p.get<double>("MinOriginalHitFraction", 0.8);
    if (maxCandidates_ == 0 || maxCandidates_ > 4 || !std::isfinite(minCoverage_) ||
        minCoverage_ < 0. || minCoverage_ > 1.)
      throw std::invalid_argument("Invalid recovery search limits");
    matcher_.DebugMode(clusmtool::CMManagerBase::kNone);
    matcher_.AddMatchAlgo(art::make_tool<clusmtool::CFloatAlgoBase>(p.get<fhicl::ParameterSet>("MatchTool")));
    matcher_.AddPriorityAlgo(art::make_tool<clusmtool::CPriorityAlgoBase>(p.get<fhicl::ParameterSet>("PriorityTool")));
  }

  RecoveryResult run(art::Event& e, const std::string& clusterLabel, const std::string& vertexLabel,
                     const HitList& originalHits, const std::vector<double>& calibration,
                     ShrRecoManager& reco, const util::GeometryUtilities& geometry) {
    RecoveryResult none;
    const auto vertices = e.getValidHandle<std::vector<recob::Vertex>>(vertexLabel);
    if (vertices->size() != 1 || originalHits.empty()) return none;
    const auto clusters = e.getValidHandle<std::vector<recob::Cluster>>(clusterLabel);
    art::FindManyP<recob::Hit> associations(clusters, e, clusterLabel);
    const auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    cluster::ClusterMaker maker;
    if (!maker.loadVertex(clock, vertices)) return none;
    const auto& wires = art::ServiceHandle<geo::WireReadout>()->Get();
    // Match the existing ClusterMaker and ProtoShowerCMTool job-level conversions.
    const auto jobClock = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    const auto properties = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(jobClock);
    const double wire2cm = wires.Plane(geo::PlaneID{0,0,0}).WirePitch();
    const double time2cm = sampling_rate(jobClock)/1000. * properties.DriftVelocity(properties.Efield(),properties.Temperature());
    const auto vertex = vertices->at(0).position();
    struct Candidate { std::size_t cluster; Split split; };
    std::array<std::vector<Candidate>, 3> candidates;
    std::vector<HitList> allHits;
    std::set<art::Ptr<recob::Hit>> uniqueHits;
    for (std::size_t c = 0; c < clusters->size(); ++c) {
      allHits.push_back(associations.at(c));
      const auto& hits = allHits.back();
      if (hits.empty()) continue;
      const unsigned int plane = hits.front()->WireID().Plane;
      if (plane >= 3) return none;
      std::vector<Hit2D> points;
      for (const auto& h : hits) {
        // Do not create ambiguous/doubled hit ownership from malformed input.
        if (h->WireID().Plane != plane || !uniqueHits.insert(h).second) return none;
        const double vw = wires.Plane(geo::PlaneID(0,0,plane)).WireCoordinate(vertex)*wire2cm + 0.15;
        const double vt = vertex.X() + trigger_offset(clock)*time2cm + plane*0.3;
        points.push_back({h->WireID().Wire*wire2cm-vw, h->PeakTime()*time2cm-vt, h->Integral()});
      }
      auto split = splitAtVertex(points, cfg_);
      if (split.valid()) candidates[plane].push_back({c, std::move(split)});
    }
    for (auto& plane : candidates) {
      std::stable_sort(plane.begin(), plane.end(), [](const Candidate& a, const Candidate& b) {
        return a.split.cutFraction < b.split.cutFraction;
      });
      if (plane.size() > maxCandidates_) plane.resize(maxCandidates_);
    }
    const std::set<art::Ptr<recob::Hit>> original(originalHits.begin(), originalHits.end());
    // Each plane can retain its existing clusters or split one candidate.
    // At most (MaxCandidatesPerPlane+1)^3-1 trials; first successful trial wins.
    for (std::size_t a = 0; a <= candidates[0].size(); ++a)
      for (std::size_t b = 0; b <= candidates[1].size(); ++b)
        for (std::size_t c = 0; c <= candidates[2].size(); ++c) {
          if (a+b+c == 0) continue;
          const std::array<std::size_t,3> choices{{a,b,c}};
          std::map<std::size_t, const Split*> selected;
          for (std::size_t pl = 0; pl < 3; ++pl) if (choices[pl]) {
            const auto& candidate = candidates[pl][choices[pl]-1];
            selected[candidate.cluster] = &candidate.split;
          }
          RecoveryResult trial;
          for (std::size_t ci = 0; ci < allHits.size(); ++ci) {
            if (allHits[ci].empty()) continue;
            const auto found = selected.find(ci);
            if (found == selected.end()) trial.clusterHits.push_back(allHits[ci]);
            else for (const auto& part : found->second->parts) {
              HitList hits;
              for (auto index : part) hits.push_back(allHits[ci].at(index));
              trial.clusterHits.push_back(std::move(hits));
            }
          }
          std::vector<cluster::Cluster> cmClusters;
          for (const auto& hits : trial.clusterHits) {
            cluster::Cluster cm;
            maker.MakeCluster(hits, cm);
            cmClusters.push_back(cm);
            const geo::PlaneID id(0,0,cm._plane);
            trial.clusters.emplace_back(
              cm._start_pt._w/wire2cm,0.,cm._start_pt._t/time2cm,0.,0.,cm._angle,0.,
              cm._end_pt._w/wire2cm,0.,cm._end_pt._t/time2cm,0.,0.,0.,0.,
              cm._sum_charge,0.,cm._sum_charge,0.,hits.size(),0.,0.,trial.clusters.size(),wires.Plane(id).View(),id);
          }
          matcher_.Reset();
          matcher_.SetClusters(cmClusters);
          matcher_.Process();
          const auto matches = matcher_.GetBookKeeper().GetResult();
          std::vector<protoshower::ProtoShower> protos;
          for (const auto& match : matches) {
            protoshower::ProtoShower proto;
            proto.Reset();
            proto._index = protos.size();
            proto._vertex = TVector3(vertex.X(),vertex.Y(),vertex.Z());
            proto.hasVertex(true);
            for (auto index : match) {
              const auto& cl = trial.clusters.at(index);
              cluster2d::Cluster2D twoD{};
              twoD._plane = cl.Plane().Plane;
              twoD._angle_2d = cl.StartAngle();
              twoD._opening_angle = 0.;
              twoD._start = util::PxHit(twoD._plane, cl.StartWire()*wire2cm,
                (cl.StartTick()-trigger_offset(clock))*time2cm,0.,0.,0.);
              twoD._end = util::PxHit(twoD._plane, cl.EndWire()*wire2cm,
                (cl.EndTick()-trigger_offset(clock))*time2cm,0.,0.,0.);
              for (const auto& hit : trial.clusterHits.at(index)) {
                const double energy = hit->Integral()*calibration.at(twoD._plane);
                twoD._hits.emplace_back(twoD._plane,hit->WireID().Wire*wire2cm,
                  (hit->PeakTime()-trigger_offset(clock))*time2cm,energy,energy,hit->PeakAmplitude());
              }
              proto._clusters.push_back(std::move(twoD));
            }
            proto.hasCluster2D(!proto._clusters.empty());
            protos.push_back(std::move(proto));
          }
          reco.SetProtoShowers(protos);
          std::vector<Shower_t> reconstructed;
          reco.Reconstruct(geometry, reconstructed);
          for (auto shower : reconstructed) if (usableShower(shower)) {
            trial.matches.push_back(matches.at(shower.fIndex));
            shower.fIndex = trial.showers.size();
            trial.showers.push_back(shower);
          }
          if (trial.showers.size() != 2) continue;
          std::set<art::Ptr<recob::Hit>> recovered;
          for (const auto& match : trial.matches) for (auto ci : match)
            recovered.insert(trial.clusterHits[ci].begin(),trial.clusterHits[ci].end());
          std::size_t retained = 0;
          for (const auto& hit : original) retained += recovered.count(hit);
          if (double(retained)/original.size() < minCoverage_) continue;
          // A split must actually contribute both pieces to the accepted pair.
          // Otherwise unrelated unmatched fragments could masquerade as recovery.
          bool contributes = false;
          for (const auto& entry : selected) {
            std::array<bool,2> used{{false,false}};
            for (std::size_t s = 0; s < 2; ++s) for (auto ci : trial.matches[s])
              for (const auto& hit : trial.clusterHits[ci])
                if (std::find(allHits[entry.first].begin(),allHits[entry.first].end(),hit) != allHits[entry.first].end())
                  used[s] = true;
            contributes = contributes || (used[0] && used[1]);
          }
          if (!contributes) continue;
          return trial;
        }
    return none;
  }
private:
  SplitConfig cfg_;
  unsigned int maxCandidates_;
  double minCoverage_;
  clusmtool::CMatchManager matcher_;
};

} // namespace recovery
} // namespace showerreco
#endif
