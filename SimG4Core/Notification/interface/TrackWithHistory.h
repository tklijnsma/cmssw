#ifndef SimG4Core_TrackWithHistory_H
#define SimG4Core_TrackWithHistory_H

#include "G4Track.hh"
#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4Allocator.hh"
#include <set>

class G4VProcess;
class G4TrackToParticleID;
/** The part of the information about a SimTrack that we need from
 *  a G4Track
 */

class TrackWithHistory {
public:
  /** The constructor is called at PreUserTrackingAction time, 
     *  when some of the information is not available yet.
     */
  TrackWithHistory(const G4Track *g4track);
  ~TrackWithHistory() {}

  inline void *operator new(size_t);
  inline void operator delete(void *TrackWithHistory);

  void save() { saved_ = true; }
  unsigned int trackID() const { return trackID_; }
  int particleID() const { return particleID_; }
  int parentID() const { return parentID_; }
  int genParticleID() const { return genParticleID_; }
  const math::XYZVectorD &momentum() const { return momentum_; }
  double totalEnergy() const { return totalEnergy_; }
  const math::XYZVectorD &vertexPosition() const { return vertexPosition_; }
  double globalTime() const { return globalTime_; }
  double localTime() const { return localTime_; }
  double properTime() const { return properTime_; }
  const G4VProcess *creatorProcess() const { return creatorProcess_; }
  double weight() const { return weight_; }
  void setTrackID(int i) { trackID_ = i; }
  void setParentID(int i) { parentID_ = i; }
  void setGenParticleID(int i) { genParticleID_ = i; }
  bool storeTrack() const { return storeTrack_; }
  bool saved() const { return saved_; }

  // Setting and getting flags for isPrimary
  void setIsPrimary() { isPrimary_ = true; }
  bool isPrimary() const { return isPrimary_; }

  // Boundary crossing variables
  void setCrossedBoundaryPosMom(int id, const math::XYZVectorD position, const math::XYZTLorentzVectorD momentum){
    crossedBoundary_ = true;
    idAtBoundary_ = id;
    positionAtBoundary_ = position;
    momentumAtBoundary_ = momentum;
    correctedMomentumAtBoundary_ = momentum;
    }
  bool crossedBoundary() const { return crossedBoundary_; }
  math::XYZVectorD getPositionAtBoundary() const {
    assertCrossedBoundary();
    return positionAtBoundary_;
    }
  math::XYZTLorentzVectorD getMomentumAtBoundary() const {
    assertCrossedBoundary();
    return momentumAtBoundary_;
    }
  math::XYZTLorentzVectorD getCorrectedMomentumAtBoundary() const {
    assertCrossedBoundary();
    return correctedMomentumAtBoundary_;
    }
  int getIDAtBoundary() const {
    assertCrossedBoundary();
    return idAtBoundary_;
    }
  void applyCorrectionToMomentumAtBoundary(unsigned int secondaryTrackID, math::XYZTLorentzVectorD secondaryMomentum){
    assertCrossedBoundary();

    if (isCorrectedByTheseTrackIDs_.count(secondaryTrackID) > 0){
      edm::LogVerbatim("DoFineCalo")
        << "Track " << trackID_ << "is already corrected for secondary " << secondaryTrackID;
      }
    else {
      math::XYZTLorentzVectorD correctedMomentumAtBoundaryBefore = correctedMomentumAtBoundary_; // Save copy for logging purposes only
      correctedMomentumAtBoundary_ -= secondaryMomentum;
      edm::LogVerbatim("DoFineCalo")
        << "Correcting " << trackID_
        << ": ("
        << correctedMomentumAtBoundaryBefore.x() << ","
        << correctedMomentumAtBoundaryBefore.y() << ","
        << correctedMomentumAtBoundaryBefore.z() << ","
        << correctedMomentumAtBoundaryBefore.e() << ")"
        << "- ("
        << secondaryMomentum.x() << ","
        << secondaryMomentum.y() << ","
        << secondaryMomentum.z() << ","
        << secondaryMomentum.e() << ")"
        << " = ("
        << correctedMomentumAtBoundary_.x() << ","
        << correctedMomentumAtBoundary_.y() << ","
        << correctedMomentumAtBoundary_.z() << ","
        << correctedMomentumAtBoundary_.e() << ")"
        ;
      isCorrectedByTheseTrackIDs_.insert(secondaryTrackID);
      }
    }
    void applyCorrectionToMomentumAtBoundary(const G4Track* secondaryTrack){
      math::XYZTLorentzVectorD secondaryMomentum(
        secondaryTrack->GetMomentum().x() / CLHEP::GeV,
        secondaryTrack->GetMomentum().y() / CLHEP::GeV,
        secondaryTrack->GetMomentum().z() / CLHEP::GeV,
        secondaryTrack->GetKineticEnergy() / CLHEP::GeV
        );
      applyCorrectionToMomentumAtBoundary(secondaryTrack->GetTrackID(), secondaryMomentum);
      }
    void applyCorrectionToMomentumAtBoundary(const TrackWithHistory* secondaryTrack){
      math::XYZTLorentzVectorD secondaryMomentum(
        secondaryTrack->momentum().x() / CLHEP::GeV,
        secondaryTrack->momentum().y() / CLHEP::GeV,
        secondaryTrack->momentum().z() / CLHEP::GeV,
        secondaryTrack->totalEnergy() / CLHEP::GeV
        );
      applyCorrectionToMomentumAtBoundary(secondaryTrack->trackID(), secondaryMomentum);
      }

//   // Getter/setter for corrected momentum at boundary. Returns ordinary momentum at boundary if not specified.
//   bool hasCorrectedMomentumAtBoundary() const {return hasCorrectedMomentumAtBoundary_;}
//   math::XYZTLorentzVectorD getCorrectedMomentumAtBoundary() const {
//     return (hasCorrectedMomentumAtBoundary_) ? correctedMomentumAtBoundary_ : getMomentumAtBoundary();
//     }
//   void setCorrectedMomentumAtBoundary(math::XYZTLorentzVectorD corrMom){
//     // Don't overwrite if new 4-mom has a higher energy
//     if (hasCorrectedMomentumAtBoundary_ && corrMom.E() > correctedMomentumAtBoundary_.E()){
// #ifdef EDM_ML_DEBUG
//       edm::LogVerbatim("DoFineCalo") << "Not overwriting correctedMomentumAtBoundary for track " << trackID_;
// #endif
//       return;
//       }
//     hasCorrectedMomentumAtBoundary_ = true;
//     correctedMomentumAtBoundary_ = corrMom;
//     edm::LogVerbatim("DoFineCalo")
//       << "Setting track " << trackID_ << " correctedMomentumAtBoundary[GeV]=("
//       << corrMom.Px() << ","
//       << corrMom.Py() << ","
//       << corrMom.Pz() << ","
//       << corrMom.E() << ")"
//       ;
//     }

  void setParentMomentumAtCreation(math::XYZTLorentzVectorD fparentMomentum) {
    hasParentMomentumAtCreation_ = true;
    parentMomentumAtCreation_ = fparentMomentum;
    }
  math::XYZTLorentzVectorD parentMomentumAtCreation() const {
    if (!hasParentMomentumAtCreation_) throw cms::Exception("Unknown", "TrackWithHistory")
        << "parentMomentumAtCreation called for track " << trackID_ << ", but it is not set";
    return parentMomentumAtCreation_;
    }

  bool passesCaloSplittingCriterion() const { return flagCaloSplittingCriterion_; }
  void setPassesCaloSplittingCriterion() { flagCaloSplittingCriterion_ = true; }


  /** Internal consistency check (optional).
     *  Method called at PostUserTrackingAction time, to check
     *  if the information is consistent with that provided
     *  to the constructor.
     */
  void checkAtEnd(const G4Track *);

private:
  unsigned int trackID_;
  int particleID_;
  int parentID_;
  int genParticleID_;
  math::XYZVectorD momentum_;
  double totalEnergy_;
  math::XYZVectorD vertexPosition_;
  double globalTime_;
  double localTime_;
  double properTime_;
  const G4VProcess *creatorProcess_;
  double weight_;
  bool storeTrack_;
  bool saved_;

  bool flagCaloSplittingCriterion_;
  bool isPrimary_;
  bool crossedBoundary_;
  int idAtBoundary_;
  math::XYZVectorD positionAtBoundary_;
  math::XYZTLorentzVectorD momentumAtBoundary_;
  // bool hasCorrectedMomentumAtBoundary_;
  math::XYZTLorentzVectorD correctedMomentumAtBoundary_;
  bool hasParentMomentumAtCreation_;
  math::XYZTLorentzVectorD parentMomentumAtCreation_;
  std::set<unsigned int> isCorrectedByTheseTrackIDs_;

  int extractGenID(const G4Track *gt) const;

  void assertCrossedBoundary() const {
    if (!crossedBoundary_){
      throw cms::Exception("Unknown", "TrackWithHistory")
        << "Assert crossed boundary failed for track " << trackID_;
      }
    }
};

extern G4ThreadLocal G4Allocator<TrackWithHistory> *fpTrackWithHistoryAllocator;

inline void *TrackWithHistory::operator new(size_t) {
  if (!fpTrackWithHistoryAllocator)
    fpTrackWithHistoryAllocator = new G4Allocator<TrackWithHistory>;
  return (void *)fpTrackWithHistoryAllocator->MallocSingle();
}

inline void TrackWithHistory::operator delete(void *aTwH) {
  fpTrackWithHistoryAllocator->FreeSingle((TrackWithHistory *)aTwH);
}

#endif
