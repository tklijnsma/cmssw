#ifndef SimG4Core_TrackInformation_H
#define SimG4Core_TrackInformation_H

#include "FWCore/Utilities/interface/Exception.h"
#include "G4VUserTrackInformation.hh"
#include "G4Allocator.hh"
#include "G4Track.hh"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class TrackInformation : public G4VUserTrackInformation {
public:
  ~TrackInformation() override {}
  inline void *operator new(size_t);
  inline void operator delete(void *TrackInformation);

  bool storeTrack() const { return storeTrack_; }
  /// can only be set to true, cannot be reset to false!
  void storeTrack(bool v) {
    if (v)
      storeTrack_ = v;
    if (v == true)
      putInHistory();
  }

  bool isPrimary() const { return isPrimary_; }
  void isPrimary(bool v) { isPrimary_ = v; }

  bool hasHits() const { return hasHits_; }
  void hasHits(bool v) { hasHits_ = v; }

  bool isGeneratedSecondary() const { return isGeneratedSecondary_; }
  void isGeneratedSecondary(bool v) { isGeneratedSecondary_ = v; }

  bool isInHistory() const { return isInHistory_; }
  void putInHistory() { isInHistory_ = true; }

  bool isAncestor() const { return flagAncestor_; }
  void setAncestor() { flagAncestor_ = true; }

  // Calo section
  int getIDonCaloSurface() const { return idOnCaloSurface_; }
  void setIDonCaloSurface(int id, int ical, int last, int pdgID, double p) {
    idOnCaloSurface_ = id;
    idCaloVolume_ = ical;
    idLastVolume_ = last;
    caloSurfaceParticlePID_ = pdgID;
    caloSurfaceParticleP_ = p;
  }
  int getIDCaloVolume() const { return idCaloVolume_; }
  int getIDLastVolume() const { return idLastVolume_; }
  bool caloIDChecked() const { return caloIDChecked_; }
  void setCaloIDChecked(bool f) { caloIDChecked_ = f; }
  int caloSurfaceParticlePID() const { return caloSurfaceParticlePID_; }
  void setCaloSurfaceParticlePID(int id) { caloSurfaceParticlePID_ = id; }
  double caloSurfaceParticleP() const { return caloSurfaceParticleP_; }
  void setCaloSurfaceParticleP(double p) { caloSurfaceParticleP_ = p; }
  int getIDfineCalo() const { return ((idFineCalo_ > 0) ? idFineCalo_ : idOnCaloSurface_); }
  void setIDfineCalo(int id) { idFineCalo_ = id; }

  bool passesCaloSplittingCriterion() const { return flagCaloSplittingCriterion_; }
  void setPassesCaloSplittingCriterion() { flagCaloSplittingCriterion_ = true; }


  // Boundary crossing variables
  void setCrossedBoundary(const G4Track* track){
    crossedBoundary_ = true;
    // Double-check units! Any conversions necessary? Is x,y,z,E fine for XYZTLorentzVectorD?
    idAtBoundary_ = track->GetTrackID();
    positionAtBoundary_ = math::XYZVectorD(
      track->GetPosition().x(),
      track->GetPosition().y(),
      track->GetPosition().z()
      );
    momentumAtBoundary_ = math::XYZTLorentzVectorD(
      track->GetMomentum().x(), track->GetMomentum().y(), track->GetMomentum().z(), track->GetKineticEnergy()
      );
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
  int getIDAtBoundary() const {
    assertCrossedBoundary();
    return idAtBoundary_;
    }

  // Generator information
  int genParticlePID() const { return genParticlePID_; }
  void setGenParticlePID(int id) { genParticlePID_ = id; }
  double genParticleP() const { return genParticleP_; }
  void setGenParticleP(double p) { genParticleP_ = p; }

  // remember the PID of particle entering the CASTOR detector. This is needed
  // in order to scale the hadronic response
  bool hasCastorHit() const { return hasCastorHit_; }
  void setCastorHitPID(const int pid) {
    hasCastorHit_ = true;
    castorHitPID_ = pid;
  }
  int getCastorHitPID() const { return castorHitPID_; }

  void Print() const override;

  void insertMomentumAtCreationSecondary(const G4Track* secondary, const G4Track* mother){
    if (momentumAtCreationSecondaryMap_.count(secondary) > 0) { return; } // If already in map, don't add again
    edm::LogVerbatim("DoFineCalo")
      << "Mother " << mother->GetTrackID()
      << ": Inserting momentumAtCreation for secondary " << secondary->GetTrackID()
      << " (address " << secondary << ")"
      << " momentumAtCreation=("
      << mother->GetMomentum().x() << ","
      << mother->GetMomentum().y() << ","
      << mother->GetMomentum().z() << ","
      << mother->GetKineticEnergy() << ")"
      << " Esecondary=" << secondary->GetKineticEnergy()
      ;
    momentumAtCreationSecondaryMap_.insert(
      std::pair<const G4Track*, math::XYZTLorentzVectorD>(
        secondary, math::XYZTLorentzVectorD(
          mother->GetMomentum().x(), mother->GetMomentum().y(), mother->GetMomentum().z(),
          mother->GetKineticEnergy()
          )
        )
      );
    }

  math::XYZTLorentzVectorD momentumAtCreation(const G4Track* secondary) const {
    auto momentumAtCreationPair = momentumAtCreationSecondaryMap_.find(secondary);
    if ( momentumAtCreationPair == momentumAtCreationSecondaryMap_.end() ) {
      throw cms::Exception("Unknown", "TrackInformation")
        << "Requested momentumAtCreation of track " << secondary->GetTrackID()
        << " with address " << secondary
        << ", but it is not a daughter of " << getIDfineCalo()
        ;
      }
    return momentumAtCreationPair->second;
    }

  math::XYZTLorentzVectorD parentMomentumAtCreation() const {
    if (!hasParentMomentumAtCreation_)
      throw cms::Exception("Unknown", "TrackInformation")
      << "Attempted to get parentMomentumAtCreation, but it is not set";
    return parentMomentumAtCreation_;
    }
  void setParentMomentumAtCreation(math::XYZTLorentzVectorD fparentMomentum){
    hasParentMomentumAtCreation_ = true;
    parentMomentumAtCreation_ = fparentMomentum;
    }

private:
  bool storeTrack_;
  bool isPrimary_;
  bool hasHits_;
  bool isGeneratedSecondary_;
  bool isInHistory_;
  bool flagAncestor_;
  int idOnCaloSurface_;
  int idCaloVolume_;
  int idLastVolume_;
  bool caloIDChecked_;
  int idFineCalo_;
  bool crossedBoundary_;
  bool idAtBoundary_;
  math::XYZVectorD positionAtBoundary_;
  math::XYZTLorentzVectorD momentumAtBoundary_;
  bool flagCaloSplittingCriterion_;

  int genParticlePID_, caloSurfaceParticlePID_;
  double genParticleP_, caloSurfaceParticleP_;

  bool hasCastorHit_;
  int castorHitPID_;

  bool hasParentMomentumAtCreation_;
  math::XYZTLorentzVectorD parentMomentumAtCreation_;
  std::map< const G4Track*, math::XYZTLorentzVectorD > momentumAtCreationSecondaryMap_;

  void assertCrossedBoundary() const {
    if (!crossedBoundary_){
      throw cms::Exception("Unknown", "TrackInformation")
        << "Assert crossed boundary failed for track "
        << getIDonCaloSurface() << " (fine: " << getIDfineCalo() << ")"
        ;
      }
    }

  // Restrict construction to friends
  TrackInformation()
      : G4VUserTrackInformation(),
        storeTrack_(false),
        isPrimary_(false),
        hasHits_(false),
        isGeneratedSecondary_(false),
        isInHistory_(false),
        flagAncestor_(false),
        idOnCaloSurface_(0),
        idCaloVolume_(-1),
        idLastVolume_(-1),
        caloIDChecked_(false),
        idFineCalo_(-1),
        crossedBoundary_(false),
        flagCaloSplittingCriterion_(false),
        genParticlePID_(-1),
        caloSurfaceParticlePID_(0),
        genParticleP_(0),
        caloSurfaceParticleP_(0),
        hasCastorHit_(false),
        castorHitPID_(0),
        hasParentMomentumAtCreation_(false)
        {}
  friend class NewTrackAction;
};

extern G4ThreadLocal G4Allocator<TrackInformation> *fpTrackInformationAllocator;

inline void *TrackInformation::operator new(size_t) {
  if (!fpTrackInformationAllocator)
    fpTrackInformationAllocator = new G4Allocator<TrackInformation>;
  return (void *)fpTrackInformationAllocator->MallocSingle();
}

inline void TrackInformation::operator delete(void *trkInfo) {
  fpTrackInformationAllocator->FreeSingle((TrackInformation *)trkInfo);
}

#endif
