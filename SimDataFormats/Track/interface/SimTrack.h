#ifndef SimTrack_H
#define SimTrack_H

#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
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
  void setCrossedBoundaryPosMom(int id, const math::XYZVectorD position, const math::XYZTLorentzVectorD momentum){
    crossedBoundary_ = true;
    idAtBoundary_ = id;
    positionAtBoundary_ = position;
    momentumAtBoundary_ = momentum;
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

private:
  int ivert;
  int igenpart;

  math::XYZVectorD tkposition;
  math::XYZTLorentzVectorD tkmomentum;

  bool crossedBoundary_;
  int idAtBoundary_;
  math::XYZVectorD positionAtBoundary_;
  math::XYZTLorentzVectorD momentumAtBoundary_;
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
