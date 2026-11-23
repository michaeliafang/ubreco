#include <iostream>
#ifndef UBRECO_SINGLESHOWERRECOVERY_H
#define UBRECO_SINGLESHOWERRECOVERY_H

#include "VertexSplit.h"
#include "HitOwnership.h"
#include "RecoveryEligibility.h"
#include "MatchHypotheses.h"
#include "SpacePointConsistency.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "canvas/Utilities/InputTag.h"
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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <chrono>
#include <memory>
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

// A seed may preserve the full cmerger collection or supply combined alternatives.
// Forced matches preserve the two main showers when attaching a third fragment.
struct RecoverySeed {
  std::vector<HitList> clusters;
  ClusterMatches forcedMatches;
  std::vector<HitList> coverageGroups;
  HitList coverageReference, detachedHits;
  double selectedHitCoverage=0.; // zero preserves the previous coverage policy
  bool mergeOnly=false;
  bool preferExisting=false;
  unsigned trialLimit=0; // zero uses the remaining shared event budget
  std::string mode="combine_split";
};

// Mirrors the pre-publication validity checks in ShrReco3D::SaveShower.
inline bool usableShower(const Shower_t& s) {
  return s.fPassedReconstruction && s.fXYZStart.Mag2() != 0. && s.fDCosStart.Mag2() != 0.;
}

class SingleShowerRecovery {
public:
  explicit SingleShowerRecovery(const fhicl::ParameterSet& p) : matcher_(3) {
    cfg_.minHits = p.get<unsigned int>("MinHits", 10);
    cfg_.minUpstreamHits = p.get<unsigned int>("MinUpstreamHits", 6);
    cfg_.maxHits = p.get<unsigned int>("MaxHits", 6000);
    cfg_.neighborRadius = p.get<double>("NeighborRadius", 0.8);
    cfg_.innerRadius = p.get<double>("InnerRadius", 0.);
    cfg_.upstreamRadius = p.get<double>("UpstreamRadius", 10.);
    cfg_.angleStepDegrees = p.get<double>("AngleStepDegrees", 2.);
    cfg_.minOpeningDegrees = p.get<double>("MinOpeningDegrees", 5.);
    cfg_.minChargeFraction = p.get<double>("MinChargeFraction", 0.20);
    cfg_.maxCutFraction = p.get<double>("MaxCutFraction", 0.15);
    cfg_.maxStartRadius = p.get<double>("MaxStartRadius", 5.);
    cfg_.minStartHits = p.get<unsigned int>("MinStartHits", 3);
    cfg_.minTrackSpan = p.get<double>("MinTrackSpan", 3.);
    cfg_.maxTrackGap = p.get<double>("MaxTrackGap", 3.0);
    cfg_.trackWidth = p.get<double>("TrackWidth", 0.6);
    cfg_.minTrackFraction = p.get<double>("MinTrackFraction", 0.0);
    cfg_.maxEdges = p.get<unsigned int>("MaxEdges", 200000);
    cfg_.validate();
    useSpacePoints_ = p.get<bool>("UseSpacePoints",false);
    if (useSpacePoints_) {
      spacePointTag_ = p.get<art::InputTag>("SpacePointProducer");
      spacePointHitTag_ = p.get<art::InputTag>("SpacePointHitProducer");
    }
    minPointWeight_ = p.get<double>("MinSpacePointWeight",6.);
    maxPointConflict_ = p.get<double>("MaxSpacePointConflict",0.25);
    if (!std::isfinite(minPointWeight_) || minPointWeight_<=0. ||
        !std::isfinite(maxPointConflict_) || maxPointConflict_<0. || maxPointConflict_>1.)
      throw std::invalid_argument("Invalid spacepoint consistency configuration");
    searchSubstantialPairs_ = p.get<bool>("SearchSubstantialClusterPairs",false);
    maxPairClusters_ = p.get<unsigned>("MaxPairClustersPerPlane",3);
    minPairChargeCoverage_ = p.get<double>("MinPairPlaneChargeCoverage",0.8);
    if (maxPairClusters_<2 || maxPairClusters_>4 || !std::isfinite(minPairChargeCoverage_) ||
        minPairChargeCoverage_<=0. || minPairChargeCoverage_>1.)
      throw std::invalid_argument("Invalid substantial cluster search configuration");
    allowTwoPlaneFallback_ = p.get<bool>("AllowTwoPlaneFallback",false);
    minTwoPlaneWeight_ = p.get<double>("MinTwoPlaneShowerWeight",3.);
    if (!std::isfinite(minTwoPlaneWeight_) || minTwoPlaneWeight_<=0.)
      throw std::invalid_argument("Invalid two-plane support threshold");
    tryAlternativePairings_ = p.get<bool>("TryAlternativePairings",false);
    pairingTool_ = art::make_tool<clusmtool::CFloatAlgoBase>(p.get<fhicl::ParameterSet>("MatchTool"));
    maxMatchTrials_ = p.get<unsigned int>("MaxMatchTrials", 4);
    if (maxMatchTrials_ == 0 || maxMatchTrials_ > 124)
      throw std::invalid_argument("MaxMatchTrials must be between 1 and 124");
    maxRecoShowers_ = p.get<unsigned int>("MaxRecoShowers", 4);
    if (maxRecoShowers_ < 2 || maxRecoShowers_ > 100)
      throw std::invalid_argument("MaxRecoShowers must be between 2 and 100");
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
                     ShrRecoManager& reco, const util::GeometryUtilities& geometry,
                     RecoveryResult* diagnostic = nullptr, bool diagnosticOnly = false,
                     double maxOutputStartDistance = 0., const RecoverySeed* seed = nullptr,
                     RecoveryBudget* sharedBudget = nullptr) {
    RecoveryBudget localBudget;
    auto& budget=sharedBudget ? *sharedBudget : localBudget;
    const bool mergeOnly=seed && seed->mergeOnly;
    const auto started = std::chrono::steady_clock::now();
    std::size_t trials = 0, matched = 0, recoShowers = 0, usablePairs = 0, coverageRejected = 0, contributionRejected = 0, budgetSkipped = 0;
    std::size_t duplicatesRemoved = 0, overlappingInputHits = 0, ownershipRejected = 0;
    const auto report = [&](const char* status) {
      const double ms = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-started).count();
      std::cout << "[SingleShowerRecovery] " << "event " << e.id() << " status=" << status
        << " trials=" << trials << " matched_pairs=" << matched
        << " reco_showers=" << recoShowers << " usable_pairs=" << usablePairs
        << " coverage_rejected=" << coverageRejected << " contribution_rejected=" << contributionRejected
        << " duplicates_removed=" << duplicatesRemoved << " overlapping_input_hits=" << overlappingInputHits
        << " ownership_rejected=" << ownershipRejected
        << " event_match_trials=" << budget.trials << " event_reco_showers=" << budget.showers
        << " budget_skipped=" << budgetSkipped << " recovery_ms=" << ms << std::endl;
    };
    RecoveryResult none;
    const auto vertices = e.getValidHandle<std::vector<recob::Vertex>>(vertexLabel);
    if (vertices->size() != 1 || (!diagnosticOnly && originalHits.empty())) { report("missing_vertex_or_hits"); return none; }
    std::vector<HitList> pointHits;
    if (useSpacePoints_ && !diagnosticOnly) {
      const auto points = e.getValidHandle<std::vector<recob::SpacePoint>>(spacePointTag_);
      const art::FindManyP<recob::Hit> links(points,e,spacePointHitTag_);
      for (std::size_t i=0; i<points->size(); ++i) {
        const auto hits = links.at(i);
        if (hits.size() >= 2) pointHits.push_back(hits);
      }
    }
    const auto clusters = e.getValidHandle<std::vector<recob::Cluster>>(clusterLabel);
    art::FindManyP<recob::Hit> associations(clusters, e, clusterLabel);
    const auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    cluster::ClusterMaker maker;
    if (!maker.loadVertex(clock, vertices)) { report("vertex_load_failed"); return none; }
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
    for (std::size_t c = 0; c < (seed ? seed->clusters.size() : clusters->size()); ++c) {
      const auto inputHits = seed ? seed->clusters.at(c) : associations.at(c);
      allHits.push_back(uniqueHitsInOrder(inputHits));
      duplicatesRemoved += inputHits.size()-allHits.back().size();
      const auto& hits = allHits.back();
      if (hits.empty()) continue;
      const unsigned int plane = hits.front()->WireID().Plane;
      if (plane >= 3) { report("invalid_plane"); return none; }
      std::vector<Hit2D> points;
      for (const auto& h : hits) {
        // Overlapping input clusters remain alternatives until output ownership is checked.
        overlappingInputHits += !uniqueHits.insert(h).second;
        if (h->WireID().Plane != plane) {
          std::cout << "[SingleShowerRecovery] event " << e.id()
            << " cluster=" << c << " plane=" << plane << " hit_key=" << h.key()
            << " hit_product=" << h.id() << " hit_plane=" << h->WireID().Plane << std::endl;
          report("mixed_plane_cluster");
          return none;
        }
        const double vw = wires.Plane(geo::PlaneID(0,0,plane)).WireCoordinate(vertex)*wire2cm + 0.15;
        const double vt = vertex.X() + trigger_offset(clock)*time2cm + plane*0.3;
        points.push_back({h->WireID().Wire*wire2cm-vw, h->PeakTime()*time2cm-vt, h->Integral()});
      }
      if (mergeOnly) continue;
      SplitDiagnostics diagnostics;
      auto split = splitAtVertex(points, cfg_, &diagnostics);
      if (points.size() >= 2*cfg_.minHits)
        std::cout << "[SingleShowerRecovery] " << "event " << e.id() << " plane=" << plane << " cluster=" << c
          << " hits=" << points.size() << " track_pairs=" << diagnostics.trackPairs
          << " balanced_pairs=" << diagnostics.balancedPairs << " connection_pass=" << diagnostics.connectedPairs
          << " valley_score=" << split.valleyScore << " valley_bands=" << split.valleyBands
          << " edge_limit=" << diagnostics.edgeLimit << " split=" << split.valid() << std::endl;
      if (split.valid()) candidates[plane].push_back({c, std::move(split)});
    }
    for (auto& plane : candidates) {
      std::stable_sort(plane.begin(), plane.end(), [](const Candidate& a, const Candidate& b) {
        return betterSplit(a.split,b.split);
      });
      if (plane.size() > maxCandidates_) plane.resize(maxCandidates_);
    }
    // Use exactly the same cluster conversion for diagnostics and reconstruction.
    const auto appendCluster = [&](RecoveryResult& output, const HitList& hits) {
      cluster::Cluster cm;
      maker.MakeCluster(hits, cm);
      const geo::PlaneID id(0,0,cm._plane);
      output.clusters.emplace_back(
        cm._start_pt._w/wire2cm,0.,cm._start_pt._t/time2cm,0.,0.,cm._angle,0.,
        cm._end_pt._w/wire2cm,0.,cm._end_pt._t/time2cm,0.,0.,0.,0.,
        cm._sum_charge,0.,cm._sum_charge,0.,hits.size(),0.,0.,output.clusters.size(),wires.Plane(id).View(),id);
      return cm;
    };
    if (diagnostic) {
      if (!seed) *diagnostic = RecoveryResult{};
      for (std::size_t plane = 0; plane < candidates.size(); ++plane)
        for (const auto& candidate : candidates[plane]) {
          const auto first = diagnostic->clusters.size();
          for (const auto& part : candidate.split.parts) {
            HitList hits;
            for (auto index : part) hits.push_back(allHits[candidate.cluster].at(index));
            appendCluster(*diagnostic,hits);
            diagnostic->clusterHits.push_back(std::move(hits));
          }
          std::cout << "[SingleShowerRecovery] event " << e.id()
            << " candidate_plane=" << plane << " source_cluster=" << candidate.cluster
            << " diagnostic_clusters=" << first << "," << first+1 << std::endl;
        }
    }
    if (diagnosticOnly) { report("diagnostic_only"); return none; }
    const std::set<art::Ptr<recob::Hit>> original(originalHits.begin(), originalHits.end());
    // Each plane can retain its existing clusters or split one candidate.
    // At most (MaxCandidatesPerPlane+1)^3-1 trials; first successful trial wins.
    std::vector<std::array<std::size_t,3>> combinations;
    for (std::size_t a = 0; a <= candidates[0].size(); ++a)
      for (std::size_t b = 0; b <= candidates[1].size(); ++b)
        for (std::size_t c = 0; c <= candidates[2].size(); ++c)
          if (a+b+c) combinations.push_back({{a,b,c}});
    // Seed recovery first tests one promising split against existing clusters.
    // Combined-all inputs still favor simultaneous multi-plane splits.
    std::stable_sort(combinations.begin(),combinations.end(),[&](const auto& a, const auto& b) {
      const auto na=std::count_if(a.begin(),a.end(),[](auto x){return x!=0;});
      const auto nb=std::count_if(b.begin(),b.end(),[](auto x){return x!=0;});
      if (!(seed && seed->preferExisting)) return na>nb;
      if (na!=nb) {
        if (na==1 || nb==1) return na==1;
        return na>nb;
      }
      if (na==1) {
        const Split *sa=nullptr, *sb=nullptr;
        for (unsigned pl=0; pl<3; ++pl) {
          if (a[pl]) sa=&candidates[pl][a[pl]-1].split;
          if (b[pl]) sb=&candidates[pl][b[pl]-1].split;
        }
        return betterSplit(*sa,*sb);
      }
      return false;
    });
    if (seed && seed->preferExisting && seed->selectedHitCoverage>0.) promoteMultiPlaneRetry(combinations);
    if (mergeOnly) combinations.push_back({{0,0,0}});
    for (const auto& choices : combinations) {
          if ((seed && seed->trialLimit && trials>=seed->trialLimit) ||
              !budget.canTry(maxMatchTrials_,maxRecoShowers_)) {
            budgetSkipped += combinations.size()-trials;
            break;
          }
          ++trials;
          ++budget.trials;
          std::map<std::size_t, const Split*> selected;
          for (std::size_t pl = 0; pl < 3; ++pl) if (choices[pl]) {
            const auto& candidate = candidates[pl][choices[pl]-1];
            selected[candidate.cluster] = &candidate.split;
          }
          RecoveryResult trial;
          std::vector<std::array<std::size_t,2>> splitChildren;
          for (std::size_t ci = 0; ci < allHits.size(); ++ci) {
            if (allHits[ci].empty()) continue;
            const auto found = selected.find(ci);
            if (found == selected.end()) trial.clusterHits.push_back(allHits[ci]);
            else {
              const auto firstChild = trial.clusterHits.size();
              splitChildren.push_back({{firstChild,firstChild+1}});
              for (const auto& part : found->second->parts) {
                HitList hits;
                for (auto index : part) hits.push_back(allHits[ci].at(index));
                trial.clusterHits.push_back(std::move(hits));
              }
            }
          }
          std::vector<cluster::Cluster> cmClusters;
          for (const auto& hits : trial.clusterHits)
            cmClusters.push_back(appendCluster(trial,hits));
          ClusterMatches initialMatches;
          if (mergeOnly) initialMatches=seed->forcedMatches;
          else {
            matcher_.Reset();
            matcher_.SetClusters(cmClusters);
            matcher_.Process();
            initialMatches=matcher_.GetBookKeeper().GetResult();
          }
          const auto& matches = initialMatches;
          std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
            << " initial_match_groups=" << matches.size() << std::endl;
          if (matches.size()<2 && !(useSpacePoints_ && searchSubstantialPairs_ && !mergeOnly)) continue;
          ++matched;
          const auto assignedFor = [&](const ClusterMatches& hypothesis) {
            std::map<art::Ptr<recob::Hit>,unsigned int> ownership;
            for (std::size_t m=0; m<hypothesis.size(); ++m)
              for (auto ci : hypothesis[m]) for (const auto& h : trial.clusterHits.at(ci))
                ownership.emplace(h,static_cast<unsigned int>(m));
            std::vector<std::vector<AssignedPointHit>> assigned;
            for (const auto& point : pointHits) {
              std::vector<AssignedPointHit> hits;
              for (const auto& h : point) {
                const auto found=ownership.find(h);
                if (found!=ownership.end()) hits.push_back({h.key(),h->WireID().Plane,found->second});
              }
              if (hits.size()>=2) assigned.push_back(std::move(hits));
            }
            return assigned;
          };
          if (useSpacePoints_) {
            const auto assigned = assignedFor(initialMatches);
            const auto detail = diagnosePoints(assigned,matches.size()==2);
            for (std::size_t k=0; k<3; ++k) {
              const char* pairs[] = {"01","02","12"};
              std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
                << " sps_pair=" << pairs[k] << " points=" << detail.pairs[k].informative
                << " weight=" << detail.pairs[k].weight << " conflict=" << detail.pairs[k].fraction() << std::endl;
            }
            std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
              << " sps_unique_points=" << detail.unique.informative << " sps_unique_weight=" << detail.unique.weight
              << " sps_unique_conflict=" << detail.unique.fraction() << std::endl;
            if (matches.size()==2) for (unsigned mask=0; mask<4; ++mask)
              std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
                << " sps_flip_mask=" << mask << " sps_flip_conflict=" << detail.flipConflict[mask] << std::endl;
          }
          std::vector<unsigned> planes;
          for (const auto& cl : trial.clusters) planes.push_back(cl.Plane().Plane);
          unsigned initialPlaneMask=0;
          for (const auto& match : initialMatches) for (auto ci : match) initialPlaneMask|=1u<<planes.at(ci);
          std::vector<MatchHypothesis> hypotheses{{initialMatches,initialPlaneMask,0}};
          if (!mergeOnly && useSpacePoints_ && tryAlternativePairings_) {
            const auto alternatives=pairingHypotheses(initialMatches,planes);
            for (const auto& h : alternatives) if (h.swapMask) hypotheses.push_back(h);
          }
          if (!mergeOnly && useSpacePoints_ && allowTwoPlaneFallback_) {
            const auto fallback=twoPlaneHypotheses(initialMatches,planes,tryAlternativePairings_);
            hypotheses.insert(hypotheses.end(),fallback.begin(),fallback.end());
          }
          if (!mergeOnly && useSpacePoints_ && searchSubstantialPairs_ && allowTwoPlaneFallback_) {
            std::vector<std::size_t> counts;
            std::vector<double> charges;
            for (const auto& hits : trial.clusterHits) {
              counts.push_back(hits.size());
              double charge=0.;
              for (const auto& hit : hits) {
                const double q=hit->Integral();
                if (!std::isfinite(q)) { charge=q; break; }
                charge+=std::max(0.,q);
              }
              charges.push_back(charge);
            }
            const auto extra=substantialPairHypotheses(planes,counts,charges,cfg_.minHits,
              maxPairClusters_,minPairChargeCoverage_);
            for (const auto& h : extra) {
              // Avoid spending reconstruction budget twice on the same pairing,
              // even if its two shower labels or cluster order are reversed.
              const auto canonical=[](ClusterMatches m) {
                for (auto& group : m) std::sort(group.begin(),group.end());
                std::sort(m.begin(),m.end()); return m;
              };
              const auto key=canonical(h.matches);
              bool duplicate=false;
              for (const auto& existing : hypotheses)
                if (canonical(existing.matches)==key) { duplicate=true; break; }
              if (!duplicate) hypotheses.push_back(h);
            }
          }
          struct ScoredHypothesis { MatchHypothesis hypothesis; PointConsistency score; };
          std::vector<ScoredHypothesis> scored;
          for (const auto& h : hypotheses) {
            if (h.matches.size()<2 || (seed && seed->preferExisting && h.matches.size()!=2)) continue;
            if (h.clusterSearch) {
              bool contributes=false;
              for (const auto& children : splitChildren) {
                const auto in=[&](unsigned m,std::size_t c) {
                  return std::find(h.matches[m].begin(),h.matches[m].end(),c)!=h.matches[m].end();
                };
                contributes=contributes || (in(0,children[0]) && in(1,children[1])) ||
                  (in(1,children[0]) && in(0,children[1]));
              }
              if (!contributes) continue;
            }
            // Alternative assignments must still satisfy the configured matching
            // tool, and no hypothesis may assign a hit to multiple showers.
            bool compatible=true;
            std::set<art::Ptr<recob::Hit>> used;
            for (const auto& match : h.matches) {
              std::set<art::Ptr<recob::Hit>> own;
              std::vector<const cluster::Cluster*> cl;
              for (auto ci : match) {
                cl.push_back(&cmClusters.at(ci));
                own.insert(trial.clusterHits.at(ci).begin(),trial.clusterHits.at(ci).end());
              }
              for (const auto& hit : own) if (!used.insert(hit).second && h.matches.size()==2) compatible=false;
              if (mergeOnly || h.swapMask || h.twoPlaneFallback) {
                pairingTool_->Reset();
                const double value=pairingTool_->Float(cl);
                if (!std::isfinite(value) || value<0.) compatible=false;
              }
            }
            if (!compatible) continue;
            // Ownership and retained-plane coverage are known before expensive
            // shower reconstruction. Do not spend the budget dropping a fragment.
            if (seed && seed->selectedHitCoverage>0.) {
              const auto planeOf=[](const auto& hit) { return hit->WireID().Plane; };
              const auto coverage=retainedPlaneCoverage(seed->coverageReference,used,h.planeMask,planeOf);
              auto augmented=used;
              augmented.insert(seed->detachedHits.begin(),seed->detachedHits.end());
              const auto withFragment=retainedPlaneCoverage(seed->coverageReference,augmented,h.planeMask,planeOf);
              if (seed->mode=="existing_cmerger" && !coverage.passes(seed->selectedHitCoverage) &&
                  withFragment.passes(seed->selectedHitCoverage)) budget.attachmentCanReachCoverage=true;
              std::cout << "[SingleShowerRecovery] event " << e.id()
                << " selected_hit_coverage=" << coverage.retained << "/" << coverage.total
                << " coverage_plane_mask=" << h.planeMask
                << " coverage_with_fragment=" << withFragment.retained << "/" << withFragment.total << std::endl;
              if (!coverage.passes(seed->selectedHitCoverage)) { ++coverageRejected; continue; }
            }
            if (seed && seed->selectedHitCoverage<=0.) {
              bool covered=true;
              for (const auto& group : seed->coverageGroups) {
                std::size_t n=0, kept=0;
                for (const auto& hit : uniqueHitsInOrder(group)) {
                  if (h.twoPlaneFallback && !(h.planeMask&(1u<<hit->WireID().Plane))) continue;
                  ++n; kept+=used.count(hit);
                }
                if (!n || double(kept)/n<minCoverage_) covered=false;
              }
              if (!covered) {
                ++coverageRejected;
                std::cout << "[SingleShowerRecovery] event " << e.id()
                  << " pre_reco_coverage_reject=1 plane_mask=" << h.planeMask << std::endl;
                continue;
              }
            }
            const auto assigned=assignedFor(h.matches);
            const auto score=pointConsistency(assigned);
            bool supported=true;
            if (mergeOnly || h.twoPlaneFallback || (seed && seed->preferExisting)) {
              const auto support=consistentShowerWeights(assigned);
              supported=useSpacePoints_ && score.weight>=minPointWeight_ && support[0]>=minTwoPlaneWeight_ && support[1]>=minTwoPlaneWeight_;
              std::array<std::array<double,3>,2> charge{};
              for (unsigned m=0; m<2; ++m) for (auto ci : h.matches[m]) {
                const auto& hits=trial.clusterHits.at(ci);
                if (hits.size()<cfg_.minHits) supported=false;
                for (const auto& hit : hits) {
                  const double q=hit->Integral();
                  if (!std::isfinite(q)) supported=false;
                  charge[m][planes[ci]]+=std::max(0.,q);
                }
              }
              for (unsigned pl=0; pl<3; ++pl) if (h.planeMask&(1u<<pl)) {
                const double total=charge[0][pl]+charge[1][pl];
                if (!std::isfinite(total) || total<=0. ||
                    std::min(charge[0][pl],charge[1][pl])<cfg_.minChargeFraction*total) supported=false;
              }
              std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
                << " fallback_plane_mask=" << h.planeMask << " fallback_swap_mask=" << h.swapMask
                << " shower_sps_weights=" << support[0] << "," << support[1]
                << " two_plane_supported=" << supported << std::endl;
            }
            const bool reject=useSpacePoints_ && score.weight>=minPointWeight_ && score.fraction()>maxPointConflict_;
            std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
              << " pairing_swap_mask=" << h.swapMask << " plane_mask=" << h.planeMask
              << " sps_weight=" << score.weight << " sps_conflict=" << score.fraction()
              << " sps_reject=" << reject << std::endl;
            // Sparse points do not justify changing the original pairing.
            if (!supported || reject || (h.swapMask && score.weight<minPointWeight_)) continue;
            scored.push_back({h,score});
          }
          std::stable_sort(scored.begin(),scored.end(),[](const auto& a,const auto& b) {
            if (a.hypothesis.twoPlaneFallback!=b.hypothesis.twoPlaneFallback)
              return !a.hypothesis.twoPlaneFallback;
            if (a.hypothesis.clusterSearch!=b.hypothesis.clusterSearch)
              return !a.hypothesis.clusterSearch;
            if (a.hypothesis.clusterSearch && b.hypothesis.clusterSearch &&
                a.hypothesis.selectedHits!=b.hypothesis.selectedHits)
              return a.hypothesis.selectedHits>b.hypothesis.selectedHits;
            return a.score.fraction()<b.score.fraction();
          });
          for (const auto& option : scored) {
            const auto& matches=option.hypothesis.matches;
            if (option.hypothesis.clusterSearch) {
              std::set<unsigned> selectedClusters;
              for (const auto& match : matches) selectedClusters.insert(match.begin(),match.end());
              std::size_t unusedHits=0, unusedClusters=0;
              for (unsigned ci=0; ci<trial.clusterHits.size(); ++ci)
                if ((option.hypothesis.planeMask&(1u<<planes[ci])) && !selectedClusters.count(ci)) {
                  ++unusedClusters; unusedHits+=trial.clusterHits[ci].size();
                }
              std::cout << "[SingleShowerRecovery] event " << e.id() << " substantial_pair_search=1 selected_clusters=";
              for (auto ci : selectedClusters) std::cout << ci << ",";
              std::cout << " unused_retained_plane_clusters=" << unusedClusters
                << " unused_hit_associations=" << unusedHits << std::endl;
            }
            trial.matches.clear();
            trial.showers.clear();
            std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
              << " reconstruct_swap_mask=" << option.hypothesis.swapMask
              << " retained_plane_mask=" << option.hypothesis.planeMask
              << " two_plane_fallback=" << option.hypothesis.twoPlaneFallback << std::endl;
            // Never truncate a match set: that could hide a third usable shower.
            if (!budget.canReconstruct(matches.size(),maxRecoShowers_)) { ++budgetSkipped; continue; }
            recoShowers += matches.size();
            budget.showers += matches.size();
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
            for (const auto& shower : reconstructed) {
              std::cout << "[SingleShowerRecovery] event " << e.id() << " trial=" << trials
                << " reco_candidate=" << shower.fIndex << " passed=" << shower.fPassedReconstruction
                << " start_mag2=" << shower.fXYZStart.Mag2() << " direction_mag2=" << shower.fDCosStart.Mag2()
                << " usable=" << usableShower(shower) << " planes=";
              for (auto ci : matches.at(shower.fIndex)) std::cout << trial.clusters.at(ci).Plane().Plane << ",";
              std::cout << std::endl;
            }
            for (auto shower : reconstructed) if (usableShower(shower)) {
              trial.matches.push_back(matches.at(shower.fIndex));
              shower.fIndex = trial.showers.size();
              trial.showers.push_back(shower);
            }
            if (trial.showers.size() != 2) continue;
            ++usablePairs;
            if (maxOutputStartDistance>0.) {
              bool near = true;
              for (const auto& shower : trial.showers) {
                const double distance=(shower.fXYZStart-TVector3(vertex.X(),vertex.Y(),vertex.Z())).Mag();
                near = near && std::isfinite(distance) && distance<=maxOutputStartDistance;
              }
              if (!near) {
                std::cout << "[SingleShowerRecovery] event " << e.id() << " recovered_start_reject=1" << std::endl;
                continue;
              }
            }
            std::array<HitList,2> outputHits;
            for (std::size_t s = 0; s < 2; ++s) for (auto ci : trial.matches[s])
              outputHits[s].insert(outputHits[s].end(),trial.clusterHits[ci].begin(),trial.clusterHits[ci].end());
            if (!disjointHits(outputHits[0],outputHits[1])) { ++ownershipRejected; continue; }
            std::set<art::Ptr<recob::Hit>> recovered;
            for (const auto& match : trial.matches) for (auto ci : match)
              recovered.insert(trial.clusterHits[ci].begin(),trial.clusterHits[ci].end());
            std::size_t retained = 0, considered = 0;
            for (const auto& hit : original) {
              if (option.hypothesis.twoPlaneFallback &&
                  !(option.hypothesis.planeMask&(1u<<hit->WireID().Plane))) continue;
              ++considered;
              retained += recovered.count(hit);
            }
            std::cout << "[SingleShowerRecovery] event " << e.id() << " coverage_retained=" << retained
              << " coverage_considered=" << considered << " coverage_omitted=" << original.size()-considered << std::endl;
            if (seed && seed->selectedHitCoverage>0.) {
              unsigned retainedPlanes=0;
              for (const auto& match : trial.matches) for (auto ci : match)
                retainedPlanes|=1u<<trial.clusters.at(ci).Plane().Plane;
              const auto coverage=retainedPlaneCoverage(seed->coverageReference,recovered,retainedPlanes,
                [](const auto& hit) { return hit->WireID().Plane; });
              std::cout << "[SingleShowerRecovery] event " << e.id()
                << " final_selected_hit_coverage=" << coverage.retained << "/" << coverage.total
                << " coverage_plane_mask=" << retainedPlanes << std::endl;
              if (!coverage.passes(seed->selectedHitCoverage)) { ++coverageRejected; continue; }
            } else if (!considered || double(retained)/considered < minCoverage_) { ++coverageRejected; continue; }
            // Require coverage of every original shower separately: a small
            // fragment must not disappear behind good aggregate coverage.
            bool groupCoverage=true;
            if (seed && seed->selectedHitCoverage<=0.) for (std::size_t g=0; g<seed->coverageGroups.size(); ++g) {
              const auto unique=uniqueHitsInOrder(seed->coverageGroups[g]);
              std::size_t n=0, kept=0;
              for (const auto& hit : unique) {
                if (option.hypothesis.twoPlaneFallback &&
                    !(option.hypothesis.planeMask&(1u<<hit->WireID().Plane))) continue;
                ++n; kept+=recovered.count(hit);
              }
              std::cout << "[SingleShowerRecovery] event " << e.id() << " coverage_group=" << g
                << " retained=" << kept << " considered=" << n << std::endl;
              if (!n || double(kept)/n<minCoverage_) groupCoverage=false;
            }
            if (!groupCoverage) { ++coverageRejected; continue; }
            // A split must actually contribute both pieces to the accepted pair.
            // Otherwise unrelated unmatched fragments could masquerade as recovery.
            bool contributes = mergeOnly;
            for (const auto& children : splitChildren) {
              const auto contains = [&](std::size_t shower, std::size_t child) {
                const auto& match = trial.matches[shower];
                return std::find(match.begin(),match.end(),child) != match.end();
              };
              contributes = contributes || (contains(0,children[0]) && contains(1,children[1]))
                                        || (contains(1,children[0]) && contains(0,children[1]));
            }
            if (!contributes) { ++contributionRejected; continue; }
            report("accepted");
            return trial;
          } // ranked matching hypotheses
        }
    report("no_accepted_pair");
    return none;
  }
private:
  bool useSpacePoints_;
  bool tryAlternativePairings_;
  bool allowTwoPlaneFallback_;
  bool searchSubstantialPairs_;
  unsigned maxPairClusters_;
  double minPairChargeCoverage_;
  double minTwoPlaneWeight_;
  std::unique_ptr<clusmtool::CFloatAlgoBase> pairingTool_;
  art::InputTag spacePointTag_, spacePointHitTag_;
  double minPointWeight_, maxPointConflict_;
  SplitConfig cfg_;
  unsigned int maxCandidates_;
  unsigned int maxRecoShowers_;
  unsigned int maxMatchTrials_;
  double minCoverage_;
  clusmtool::CMatchManager matcher_;
};

} // namespace recovery
} // namespace showerreco
#endif
