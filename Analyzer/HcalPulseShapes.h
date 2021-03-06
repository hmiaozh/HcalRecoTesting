#ifndef HCALPULSESHAPES_H
#define HCALPULSESHAPES_H 1

#include <map>
#include <vector>
#include "HcalPulseShape.h"


/** \class HcalPulseShapes
  *  
  * \author J. Mans - Minnesota
  */
class HcalPulseShapes {
public:
  typedef HcalPulseShape Shape;
  HcalPulseShapes();
  ~HcalPulseShapes();
  // only needed if you'll be geting shapes by DetId
 // void beginRun(edm::EventSetup const & es);
  void endRun();

  const Shape& hbShape() const { return hpdShape_; }
  const Shape& heShape() const { return hpdShape_; }
  const Shape& hfShape() const { return hfShape_; }
  const Shape& hoShape(bool sipm=false) const { return sipm ? siPMShape_ : hpdShape_; }
  //  return Shpape for given shapeType.
  const Shape& getShape(int shapeType) const;
  /// automatically figures out which shape to return
  //const Shape& shape(const HcalDetId & detId) const;
  //const Shape& shapeForReco(const HcalDetId & detId) const;
  /// in case of conditions problems
  //const Shape& defaultShape(const HcalDetId & detId) const;
private:
  void computeHPDShape(float, float, float, float, float ,
                       float, float, float, Shape&);
  // void computeHPDShape();
  void computeHFShape();
  void computeSiPMShape();
  Shape hpdShape_, hfShape_, siPMShape_;
  Shape hpdShape_v2, hpdShapeMC_v2;
  Shape hpdShape_v3, hpdShapeMC_v3;
  Shape hpdBV30Shape_v2, hpdBV30ShapeMC_v2;
//   HcalMCParams * theMCParams;
//   const HcalTopology * theTopology;
//   HcalRecoParams * theRecoParams;
  typedef std::map<int, const Shape *> ShapeMap;
  ShapeMap theShapes;

};
#endif
