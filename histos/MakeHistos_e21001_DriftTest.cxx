
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TBank88.h"
#include "TS800.h"
#include "GCutG.h"

#include "TChannel.h"
#include "GValue.h"

#define Q1 15
#define Q2 7
#define Q3 11
#define Q4 1
#define Q5 22
#define Q6 14
#define Q7 12
#define Q8 6
#define Q9 21


std::map<int,int> HoleQMap;
std::map<int,std::string> LayerMap;

void InitMap() {
  HoleQMap[Q1] = 1;
  HoleQMap[Q2] = 2;
  HoleQMap[Q3] = 3;
  HoleQMap[Q4] = 4;
  HoleQMap[Q5] = 5;
  HoleQMap[Q6] = 6;
  HoleQMap[Q7] = 7;
  HoleQMap[Q8] = 8;
  HoleQMap[Q9] = 9;

  LayerMap[0] = "alpha";
  LayerMap[1] = "beta";
  LayerMap[2] = "gamma";
  LayerMap[3] = "delta";
  LayerMap[4] = "epsilon";
  LayerMap[5] = "phi";

}

#define INTEGRATION 128.0

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.

bool HandleGretina(TRuntimeObjects &obj) {

   TGretina *gretina = obj.GetDetector<TGretina>();
   //   TS800 *s800       = obj.GetDetector<TS800>();

   if(!gretina)
     return false;
   
   std::string dirname = "gretina";

   Int_t    energyNChannels = 2000;
   Double_t energyLlim = 0.;
   Double_t energyUlim = 8000.;
   
   for(int x=0;x<gretina->Size();x++) {
     TGretinaHit hit = gretina->GetGretinaHit(x);
     std::string histname = "energy";
     obj.FillHistogram(dirname, histname,
		       energyNChannels*4, energyLlim, energyUlim,
		       hit.GetCoreEnergy());

     histname = "hole";
     obj.FillHistogram(dirname, histname,
		       30, 0, 30,
		       hit.GetHoleNumber());

     histname = "overview";
     obj.FillHistogram(dirname, histname,
		       energyNChannels*4, energyLlim, energyUlim,
		       hit.GetCoreEnergy(),
		       120, 0, 120,
		       hit.GetCrystalId());
       
     // e21001: this crystal was shut down due to noise issues.
     if(hit.GetCrystalId()==78)
       continue;

     histname = Form("crystal_%d_energy",hit.GetCrystalId());
     obj.FillHistogram(dirname, histname,
		       energyNChannels*4, energyLlim, energyUlim,
		       hit.GetCoreEnergy());
     histname = Form("crystal_%d_energy_time",hit.GetCrystalId());
     obj.FillHistogram(dirname, histname,
		       4000,0,4000,hit.GetTime()*1e-8,
		       energyNChannels*4, energyLlim, energyUlim,
		       hit.GetCoreEnergy());

   }


  return true;
}

extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  //std::cout << "---------------------------------" <<std::endl;
  //std::cout << " At the beginning" << std::endl;
  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  //  TS800 *s800       = obj.GetDetector<TS800>();
  //std::cout << " Dets Gotten" << std::endl;
  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();

  std::string histname = "";
  std::string dirname  = "";

  //  if(s800) {
    
    if(gretina) {
      
      HandleGretina(obj);


    }
    //  }

  if(numobj!=list->GetSize())
    list->Sort();

}

