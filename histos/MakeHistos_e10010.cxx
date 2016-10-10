
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
#include "TBank29.h"
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

//#define BETA .37

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

GCutG *incoming_p43=0;
GCutG *incoming_s44=0;
GCutG *incoming_cl46=0;

GCutG *outgoing_s44=0;

GCutG *gt_time=0;

bool OutgoingS800(TRuntimeObjects &obj, TS800 *s800, GCutG *outgoing) {

   if(!s800)
     return false;
   
   std::string dirname = "S800";
   std::string histname;
   if(outgoing) {
     //     if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),s800->GetIonChamber().Charge()) )
     if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),
			    s800->GetIonChamber().GetSum()) )
       return false;
     histname = Form("incoming_%s", outgoing->GetName());
   } else {
     histname = "incoming";
   }
   obj.FillHistogram(dirname, histname,
		     1000, -2600, -1600,
		     s800->GetOBJ_E1Raw_MESY(),
		     1000, 2000, 3000,
		     s800->GetXF_E1Raw_MESY());

   return true;
}


bool IncomingS800(TRuntimeObjects &obj,TS800 *s800,GCutG *incoming) {
 
   if(!s800)
     return false;
   
   std::string dirname = "S800";
   std::string histname;
   if(incoming) {
     if( !incoming->IsInside(s800->GetOBJ_E1Raw_MESY(),
			     s800->GetXF_E1Raw_MESY()) )
       return false;
     histname = Form("outgoing_%s", incoming->GetName());
   } else {
     histname = "outgoing";
   }
   obj.FillHistogram(dirname, histname,
		     1000, -2500, -1500,
		     s800->GetCorrTOF_OBJ_MESY(),
		     1000, 30000, 70000,
		     s800->GetIonChamber().GetSum());
   //		     s800->GetIonChamber().Charge());

   // std::cout << "s800->GetCorrTOF_OBJ_MESY() = "
   // 	     << s800->GetCorrTOF_OBJ_MESY()
   // 	     << "s800->GetIonChamber().Charge() = "
   // 	     << s800->GetIonChamber().Charge()
   // 	     << std::endl;
   
   //   TMTof &mtof = s800->GetMTof();
   //histname = "RefSize";
   //obj.FillHistogram(dirname,histname,20,0,0,mtof.RefSize(),
   //                                   20,0,0,mtof.E1UpSize());
   
   //for(int i=0;i<mtof.E1UpSize();i++) {
   //  if(i<mtof.RefSize()) {
   //  //for(int j=0;j<mtof.E1UpSize();j++) {
   //    histname = "E1up";
   //    obj.FillHistogram(dirname,histname,2000,0,0,mtof.fE1Up.at(i)-mtof.fRef.at(i));
   //  }
   //  if(i<mtof.RfSize()) {
   //  //for(int j=0;j<mtof.RfSize();j++) {
   //    histname = "Rf";
   //    obj.FillHistogram(dirname,histname,5000,0,0,mtof.fE1Up.at(i)-mtof.fRf.at(i));
   //  }
   //  if(i<mtof.ObjSize()) {
   //  //for(int j=0;j<mtof.ObjSize();j++) {
   //    histname = "Obj";
   //    obj.FillHistogram(dirname,histname,5000,0,5000,mtof.fE1Up.at(i)-mtof.fObj.at(i));
   //  }
   //}

// Time corrections
//  histname = "time_x";
//  obj.FillHistogram(dirname,histname,2000,0,0,s800->GetCorrTOF_OBJ_MESY(),
//                                     2000,0,0,s800->GetCrdc(0).GetDispersiveX());
//
//  histname = "time_afp";
//  obj.FillHistogram(dirname,histname,2000,0,0,s800->GetCorrTOF_OBJ_MESY(),
//                                     2000,0,0,s800->GetAFP());
//
//  histname = "time_charge";
//  obj.FillHistogram(dirname,histname,2000,0,0,s800->GetCorrTOF_OBJ_MESY(),
//                                     2000,0,0,s800->GetIonChamber().Charge());

  return true;

}

bool DTA(TS800* s800, GCutG *incoming, GCutG *outgoing){

  if(!s800)
    return false;
  
  std::string dirname = "S800";
  std::string histname = Form("dta_%s_%s",
			      incoming->GetName(),
			      outgoing->GetName());

  if(!incoming || !outgoing)
    return false;

  if(!incoming->IsInside(s800->GetOBJ_E1Raw_MESY(),
			 s800->GetXF_E1Raw_MESY()) )
     return false;

  if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),
			 s800->GetIonChamber().GetSum()) )
    return false;
  
}

bool TriggerRegister(TRuntimeObjects &obj, GCutG *incoming, GCutG *outgoing){

  TGretina *gretina = obj.GetDetector<TGretina>();
  TS800       *s800 = obj.GetDetector<TS800>();

  if(!gretina || !s800)
    return false;
  
  std::string dirname = "S800";
  std::string histname;

  histname = "trigger_bit";
  unsigned short bits = s800->GetTrigger().GetRegistr();
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram(dirname, histname, 20, 0, 20, j);
  }
  histname = "trigger_raw";
  obj.FillHistogram(dirname, histname,
		    20, 0, 20,
		    s800->GetTrigger().GetRegistr());

  if(!incoming || !outgoing)
    return false;

  if(!incoming->IsInside(s800->GetOBJ_E1Raw_MESY(),
			 s800->GetXF_E1Raw_MESY()) )
     return false;

  if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),
			 s800->GetIonChamber().GetSum()) )
    return false;

  histname = "trigger_bit_PID";
  for(int j=0;j<16;j++) {
    if(((bits>>j)&0x0001))
      obj.FillHistogram(dirname, histname, 20, 0, 20, j);
  }
  histname = "trigger_raw_PID";
  obj.FillHistogram(dirname, histname,
		    20, 0, 20,
		    s800->GetTrigger().GetRegistr());
  return true;
  
}

bool HandleGretina(TRuntimeObjects &obj,GCutG *incoming,
		   GCutG *outgoing, GCutG *gt_time) {

   if(!incoming || !outgoing)
     return false;

   TGretina *gretina = obj.GetDetector<TGretina>();
   TS800 *s800       = obj.GetDetector<TS800>();

   if(!gretina || !s800)
     return false;

   if(!((s800->GetTrigger().GetRegistr()&0x0002)>>1))
     return false;

   if(!incoming->IsInside(s800->GetOBJ_E1Raw_MESY(),
			  s800->GetXF_E1Raw_MESY()) )
     return false;

   //   if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),
   //                          s800->GetIonChamber().Charge()) )
   if(!outgoing->IsInside(s800->GetCorrTOF_OBJ_MESY(),
			  s800->GetIonChamber().GetSum()) )
     return false;

   
   
   std::string dirname = "gretina";
   
   double beta = GValue::Value("BETA");
   
   for(int x=0;x<gretina->Size();x++) {
     TGretinaHit hit = gretina->GetGretinaHit(x);

     std::string histname = Form("doppler_%s",outgoing->GetName());
     obj.FillHistogram(dirname, histname,
		       4000, 0, 4000,
		       hit.GetDoppler(beta, 0));

     histname = Form("doppler_s800_%s",outgoing->GetName());
     TVector3 track = s800->Track();
     obj.FillHistogram(dirname, histname,
		       4000, 0, 4000,
		       hit.GetDoppler(beta, &track));

     histname = Form("doppler_theta_%s_t",outgoing->GetName());
     obj.FillHistogram(dirname, histname,
		       100, 0, TMath::Pi(),
		       hit.GetTheta(),
		       1000, 0, 4000,
		       hit.GetDoppler(beta, 0));

     // Time cut
     if( gt_time && !gt_time->IsInside(s800->GetTimestamp()-hit.GetTime(),
				       hit.GetCoreEnergy()) )
       return false;
  
     histname = Form("doppler_%s_t",outgoing->GetName());
     obj.FillHistogram(dirname, histname,
		       4000, 0, 4000,
		       hit.GetDoppler(beta, 0));

     histname = Form("doppler_s800_%s_t",outgoing->GetName());
     obj.FillHistogram(dirname, histname,
		       4000, 0, 4000,
		       hit.GetDoppler(beta, &track));

     histname = Form("doppler_theta_%s_t",outgoing->GetName());
     obj.FillHistogram(dirname, histname,
		       100, 0, TMath::Pi(),
		       hit.GetTheta(),
		       1000, 0, 4000,
		       hit.GetDoppler(beta, 0));
   }


  return true;
}


extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  //std::cout << "---------------------------------" <<std::endl;
  //std::cout << " At the beginning" << std::endl;
  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TS800 *s800       = obj.GetDetector<TS800>();
  //std::cout << " Dets Gotten" << std::endl;
  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();

  TList *gates = &(obj.GetGates());
  if(!incoming_p43) {
    incoming_p43 = (GCutG*)gates->FindObject("incoming_p43");
  }
  if(!incoming_s44) {
    incoming_s44 = (GCutG*)gates->FindObject("in44S");
  }
  if(!incoming_cl46) {
    incoming_cl46 = (GCutG*)gates->FindObject("incoming_cl46");
  }
  
  if(!outgoing_s44) {
    outgoing_s44 = (GCutG*)gates->FindObject("out44S");
  }

  if(!gt_time) {
    gt_time = (GCutG*)gates->FindObject("gt_time");
  }

  std::string histname = "";
  std::string dirname  = "";

  if(bank29) {
    for(unsigned int x=0;x<bank29->Size();x++) {
      TMode3Hit &hit = (TMode3Hit&)bank29->GetHit(x);
      dirname = "Bank29";
      histname = Form("bank29_%i",hit.GetChannel());
      obj.FillHistogram(dirname,histname,
          16000,0,64000,hit.Charge());
    }
    if(s800) {
      histname = "S800_Bank29_time";
      dirname  = "Bank29";
      obj.FillHistogram(dirname,histname,
          200,-200,200,bank29->Timestamp()-s800->Timestamp());
    }
  }

  if(s800) {
    dirname = "S800";

    IncomingS800(obj, s800, 0);
    IncomingS800(obj, s800, incoming_p43);
    IncomingS800(obj, s800, incoming_s44);
    IncomingS800(obj, s800, incoming_cl46);

    OutgoingS800(obj, s800, 0);
    OutgoingS800(obj, s800, outgoing_s44);

    DTA(obj, incoming_s44, outgoing_s44);
    
    if(gretina) {
      
      HandleGretina(obj, incoming_s44, outgoing_s44, gt_time);

      dirname = "gretina";
      for(unsigned int i=0;i<gretina->Size();i++) {
	TGretinaHit hit = gretina->GetGretinaHit(i);
	histname = "dtime_all";
	obj.FillHistogram(dirname,histname,
			  500, -250, 250,
			  s800->GetTimestamp()-hit.GetTimestamp(),
			  1000, 0, 4000,
			  hit.GetCoreEnergy());
	histname = "dtimet0_all";
	obj.FillHistogram(dirname,histname,
			  500, -250, 250,
			  s800->GetTimestamp()-hit.GetTime(),
			  1000, 0, 4000,
			  hit.GetCoreEnergy());
	unsigned short bits = s800->GetTrigger().GetRegistr();
	for(int j=0;j<5;j++) {
	  if((bits>>j)&0x0001) {
	    histname = Form("dtime_all_reg%i",j);
	    obj.FillHistogram(dirname,histname,
			      500, -250, 250,
			      s800->GetTimestamp()-hit.GetTimestamp(),
			      500, 0, 4000,
			      hit.GetCoreEnergy());
	    histname = Form("dtimet0_reg%i",j);
	    obj.FillHistogram(dirname,histname,
			      500, 0, 500,
			      s800->GetTimestamp()-hit.GetTime(),
			      500, 0, 4000,
			      hit.GetCoreEnergy());
	  }
	}
      }
    }
    
   // histname = "TDC_vs_DispX";
   // obj.FillHistogram(dirname,histname,
   //     4000,-4000,4000,s800->GetCorrTOF_OBJ(),
   //     600,-300,300,s800->GetCrdc(0).GetDispersiveX());

   // histname ="PID_TDC";
   // obj.FillHistogram(dirname,histname,
   //     2000,0,0,s800->GetCorrTOF_OBJ(),
   //     2000,0,0,s800->GetIonChamber().Charge());

    //double delta_t = s800->GetScint().GetTimeUp()-s800->GetTof().GetOBJ();
    //if(delta_t>E1_TDC_low && delta_t<E1_TDC_high){
    //  histname ="PID_TDC_Prompt";
    //  obj.FillHistogram(dirname,histname,
    //      8000,-4000,4000,s800->GetCorrTOF_OBJ(),
    //      2000,10000,60000,s800->GetIonChamber().Charge());
    //}

    //std::cout << " In gret loop + s800 before haspids" << std::endl;
    //std::cout << " In gret loop + s800 before haspids" << std::endl;
    //histname = "E1_m_TDC";
    //obj.FillHistogram(histname,8000,-8000,8000,s800->GetScint().GetTimeUp()-s800->GetTof().GetOBJ());

    //histname = "E1Raw";
    //obj.FillHistogram(histname,2000,0,8000,s800->GetScint().GetTimeUp());

    //histname = "M_E1Raw";
    //for(int z=0;z<s800->GetMTof().E1UpSize();z++)
    //  obj.FillHistogram(histname,2000,0,8000,s800->GetMTof().fE1Up[z]);

    histname = "CRDC1_X";
    dirname  = "CRDC";
    obj.FillHistogram(dirname,histname,800,-400,400,
		      s800->GetCrdc(0).GetDispersiveX());

    histname = "CRDC2_X";
    obj.FillHistogram(dirname,histname,800,-400,400,
		      s800->GetCrdc(1).GetDispersiveX());

  }
  //std::cout << " After Gret Calorimeter" << std::endl;
  if(numobj!=list->GetSize())
    list->Sort();
  //std::cout << " end" << std::endl;
}

