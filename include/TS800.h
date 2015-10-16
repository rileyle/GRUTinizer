#ifndef _TSEIGHTHUNDRAD_H_
#define _TSEIGHTHUNDRAD_H_

#include "TClonesArray.h"

#include "TDetector.h"

#include "TS800Hit.h"
//class TTrigger;

class TS800 : public TDetector {
public:
  TS800();
  virtual ~TS800();


  //////////////////////////////////////////////
  virtual void InsertHit(const TDetectorHit&);
  virtual TDetectorHit& GetHit(int i);
  virtual int Size();
  /////////////////////////////////////////////
  
  void SetEventCounter(Long_t event) { fEventCounter = event; }
  
  Long_t GetEventCounter() { return fEventCounter;}
  Long_t GetTimestamp()  { Timestamp(); }


  virtual void Copy(TObject& obj)        const;
  //virtual void Print(Option_t *opt = "") const;
  virtual void Clear(Option_t* opt = "");


private:
  virtual int  BuildHits();
  
  bool HandleTrigPacket(unsigned short*,int);     //!
  bool HandleTOFPacket(unsigned short*,int);     //!
  //bool HandleScintPacket(unsigned short*,int); //!
  //bool HandleFPStPacket(char*,int);     //!
  //bool HandleIonCPacket(unsigned short*,int);     //!
  //bool HandleCRDCPacket(char*,int);     //!
  //bool HandlePPACPacket(char*,int);     //!
  //bool HandleHODOPacket(char*,int);     //!
  //bool HandleMTDCPacket(char*,int);     //!

  //TScintillatr scint[3];
  TTrigger     trigger;
  TTof         tof;
  //TIonChamber  ion;
  //TCrdc        crdc[2];
  //THodoscope   hodo[32];
  //TMultiHitTof multi_tof;
 
  Long_t fEventCounter;

  ClassDef(TS800,1);
};

//S8OO PIN
//  -- e
//  -- t
//  -- ecorr
//  -- dE
//

//S800 PPAC
//  -- u
//  -- d
//  -- l
//  -- r
//  -- a
//  -- t
//  -- sumx
//  -- sumy
//  -- x
//  -- y
//  --
//  -- z
//





//S800 Focal Plane 
//  -- FpScintillator E1
//  -- FpScintillator E2  <-- most likely does not exist
//  -- FpScintillator E3  <-- most likely does not exist
//  -- FpIonChamber
//  -- FpCRDC  crdc1
//  -- FpCRDC  crdc2
//  -- FpTrack track
//  -- FpHodo  hodo


//S800 Target 
//  -- PPac ppac2 (why 2?) <-- class
//  -- PIN  pin1 <-- class
//  -- PIN  pin2 <-- class
//  -- E 

//S800 Image
//  -- Classic PPAC ppac1
//  -- Classic PPAC ppac2
//  -- double x
//  -- double theta 
//  -- double y
//  -- double phi
//  -- double dp
//  --        z
//  --        gap
//  --        dpphi
//  --
//  -- bool trigger
//  -- char *gate

//S800 Object
//  -- PIN   <-- class

//S800 Trigger
//  -- registr
//  -- s800
//  -- external1
//  -- external2
//  -- secondary 
//  -- timestamp

//S800 Time Of Flight
//  -- rf
//  -- obj
//  -- xfp
//  -- xfp_obj
//  -- rfe1
//  -- obje1
//  -- xfpe1
//  -- obje2
//  -- xfpe2
//  -- tac_obj
//  -- tac_obje1
//  -- tac_xfp
//  -- tac_xfpe1
//  -- tac_xfpobj
//  -- tacobj_obj
//  -- tacxfp_xfp
//  --
//  -- mesyhit
//  -- me1up
//  -- me1dwn
//  -- mobj
//  -- mxfp
//  -- mrf
//  -- mcrdc1
//  -- mcrdc2
//  -- mhoth
//  -- mref
//  -- mobje1
//  -- mxfpe1
//  -- mrfe1
//  -- 
//  -- diaor
//  -- dia1
//  -- dia2
//  -- dia3
//  -- dia4
//  -- dia1RF
//  -- dia2RF
//  -- dia3RF
//  -- dia4RF
//  -- diaRF
//  -- dia1Cor
//  -- dia2Cor
//  -- dia3Cor
//  -- dia4Cor
//  -- diaCor
//  --
//  -- obj_shift
//  -- xfp_shift
//  -- rf_shift
//  -- tac_obj_shift
//  -- tac_xfp_shift
//  -- 
//  --              // sc800.h +1253
//


//S800 LaBr
//  -- eraw
//  -- ecal 
//  -- traw
//  -- tcal
//  -- e_a0
//  -- e_a1
//  -- e_a2
//  -- t_a0
//
//
//




#endif 
