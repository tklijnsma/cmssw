#ifndef SimTrack_H
#define SimTrack_H

#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
// #include "SimG4Core/Notification/interface/G4SimTrack.h"
#include "FWCore/Utilities/interface/Exception.h"

class SimTrack : public CoreSimTrack {
public:
  typedef CoreSimTrack Core;

  /// constructor
  SimTrack();
  SimTrack(int ipart, const math::XYZTLorentzVectorD& p);

  /// full constructor (pdg type, momentum, time,
  /// index of parent vertex in final vector
  /// index of corresponding gen part in final vector)
  SimTrack(int ipart, const math::XYZTLorentzVectorD& p, int iv, int ig);

  SimTrack(int ipart,
           const math::XYZTLorentzVectorD& p,
           int iv,
           int ig,
           const math::XYZVectorD& tkp,
           const math::XYZTLorentzVectorD& tkm);

  /// constructor from transient
  SimTrack(const CoreSimTrack& t, int iv, int ig);

  /// index of the vertex in the Event container (-1 if no vertex)
  int vertIndex() const { return ivert; }
  bool noVertex() const { return ivert == -1; }

  /// index of the corresponding Generator particle in the Event container (-1 if no Genpart)
  int genpartIndex() const { return igenpart; }
  bool noGenpart() const { return igenpart == -1; }

  const math::XYZVectorD& trackerSurfacePosition() const { return tkposition; }

  const math::XYZTLorentzVectorD& trackerSurfaceMomentum() const { return tkmomentum; }

  inline void setTkPosition(const math::XYZVectorD& pos) { tkposition = pos; }

  inline void setTkMomentum(const math::XYZTLorentzVectorD& mom) { tkmomentum = mom; }

  inline void setVertexIndex(const int v) { ivert = v; }

  // Boundary crossing variables
  // void setCrossedBoundaryPosMom(int id, const math::XYZVectorD position, const math::XYZTLorentzVectorD momentum){
  //   crossedBoundary_ = true;
  //   idAtBoundary_ = id;
  //   positionAtBoundary_ = position;
  //   momentumAtBoundary_ = momentum;
  //   }
  // void copyCrossedBoundaryVars(const G4SimTrack* track){
  //   if (track->crossedBoundary()){
  //     crossedBoundary_ = track->crossedBoundary();
  //     idAtBoundary_ = track->getIDAtBoundary();
  //     positionAtBoundary_ = track->getPositionAtBoundary();
  //     momentumAtBoundary_ = track->getMomentumAtBoundary();
  //     correctedMomentumAtBoundary_ = track->getCorrectedMomentumAtBoundary();
  //     }
  //   }
  void setCrossedBoundaryVars(
    bool crossedBoundary,
    int idAtBoundary,
    math::XYZVectorD positionAtBoundary,
    math::XYZTLorentzVectorD momentumAtBoundary,
    math::XYZTLorentzVectorD correctedMomentumAtBoundary
    ){
    crossedBoundary_ = crossedBoundary;
    idAtBoundary_ = idAtBoundary;
    positionAtBoundary_ = positionAtBoundary;
    momentumAtBoundary_ = momentumAtBoundary;
    correctedMomentumAtBoundary_ = correctedMomentumAtBoundary;
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
  // // Getter/setter for corrected momentum at boundary. Returns ordinary momentum at boundary if not specified.
  // bool hasCorrectedMomentumAtBoundary() const {return hasCorrectedMomentumAtBoundary_;}
  // math::XYZTLorentzVectorD getCorrectedMomentumAtBoundary() const {
  //   return (hasCorrectedMomentumAtBoundary_) ? correctedMomentumAtBoundary_ : getMomentumAtBoundary();
  //   }
  // void setCorrectedMomentumAtBoundary(math::XYZTLorentzVectorD corrMom){
  //   hasCorrectedMomentumAtBoundary_ = true;
  //   correctedMomentumAtBoundary_ = corrMom;
  //   }

private:
  int ivert;
  int igenpart;

  math::XYZVectorD tkposition;
  math::XYZTLorentzVectorD tkmomentum;

  bool crossedBoundary_;
  int idAtBoundary_;
  math::XYZVectorD positionAtBoundary_;
  math::XYZTLorentzVectorD momentumAtBoundary_;
  // bool hasCorrectedMomentumAtBoundary_;
  math::XYZTLorentzVectorD correctedMomentumAtBoundary_;
  
  void assertCrossedBoundary() const {
    if (!crossedBoundary_){
      throw cms::Exception("Unknown", "SimTrack")
        << "Assert crossed boundary failed for track " << trackId();
      }
    }
};

#include <iosfwd>
std::ostream& operator<<(std::ostream& o, const SimTrack& t);

#endif
