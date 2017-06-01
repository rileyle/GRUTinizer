#ifndef TGRETSIMHIT_H
#define TGRETSIMHIT_H

#include <TObject.h>
#include <Rtypes.h>
#include <TVector3.h>
#include <TMath.h>

#include <cmath>

#include "TGEBEvent.h"
#include "TDetectorHit.h"

#define MAXHPGESEGMENTS 36

class TGretSimHit : public TDetectorHit {

public:
  TGretSimHit();
  ~TGretSimHit();

  void Copy(TObject& obj) const;

  //  void BuildFrom(const TRawEvent::GEBBankType1& raw); // ??

  virtual Int_t Charge()        const { return fEnergy;  }

  const char *GetName() const;



  void  Print(Option_t *opt="") const;
  void  Clear(Option_t *opt="");

  double GetEn() const   { return fEnergy; }
  double GetBeta() const { return fBeta; }
  /* double GetX()  const   { return fPosit.X(); } */
  /* double GetY()  const   { return fPosit.Y(); } */
  /* double GetZ()  const   { return fPosit.Z(); } */
  double GetX()  const   { return fInteraction.X(); }
  double GetY()  const   { return fInteraction.Y(); }
  double GetZ()  const   { return fInteraction.Z(); }
  
  double GetDoppler(const TVector3 *vec=0) {

    if(vec==0) {
      vec = &BeamUnitVec;
    }
    double tmp = 0.0;
    double gamma = 1/(sqrt(1-pow(GetBeta(),2)));
    tmp = fEnergy*gamma *(1 - GetBeta()*TMath::Cos(fPosit.Angle(*vec)));
    return tmp;
  }



  double GetPhi() const {
    double phi = fPosit.Phi();
    //    if(phi<0) {
    //      return TMath::TwoPi()+phi;
    //    } else {
    //    return phi;
    //    }
    if (phi>TMath::PiOver2())
      phi -= 3.*TMath::PiOver2();
    else
      phi += TMath::PiOver2();
    return phi;
  }
  double GetTheta() const { return fPosit.Theta(); }
  
  //void SetPosition(TVector3 &vec) { fCorePosition = vec; }

  double   fEnergy;
  TVector3 fPosit;
  TVector3 fInteraction;
  double   fBeta;


private:
  ClassDef(TGretSimHit,1)
};


#endif
