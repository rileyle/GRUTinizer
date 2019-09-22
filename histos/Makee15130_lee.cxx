#include "TRuntimeObjects.h"

#include <map>
#include <iostream>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TBank88.h"
#include "TUML.h"

#include "TChannel.h"
#include "GValue.h"
#include "GCutG.h"


//quads as of June 2019.
#define Q1 15
#define Q2 7
#define Q3 11
#define Q4 16 
#define Q5 8
#define Q6 14
#define Q7 12
#define Q8 17
#define Q9 19
#define Q10 6
#define Q11 9
//#define Q12 20

//#define BETA .37

std::map<int,int> HoleQMap;
std::map<int,std::string> LayerMap;

bool map_inited=false;

std::multimap<long, TGretinaHit> gretina_map;
std::multimap<long, TUML> uml_map;
std::multimap<long, TMode3Hit> bank88_map;

//GCutG *e1_e2=0;

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
  HoleQMap[Q10] = 10;
  HoleQMap[Q11] = 11;
  //  HoleQMap[Q12] = 12;
  LayerMap[0] = "alpha";
  LayerMap[1] = "beta";
  LayerMap[2] = "gamma";
  LayerMap[3] = "delta";
  LayerMap[4] = "epsilon";
  LayerMap[5] = "phi";
}



void HandleUML(TRuntimeObjects& obj) {

  TUML *uml = obj.GetDetector<TUML>();

  for(unsigned int i=0;i<uml->fSssd.size();i++) {
    TUMLHit hit = uml->fSssd.at(i);
    if(hit.Charge()>10) 
      obj.FillHistogram("strips",20,0,20,hit.GetChannel()-15); // 1 -16?
    obj.FillHistogram("strips_aligned",6400,0,64000,hit.GetEnergy(),
        20,0,20,hit.GetChannel()-15); // 1 -16?
  }
  
  for(unsigned int i=0;i<uml->Size();i++) {
    TUMLHit hit = uml->GetUMLHit(i);
    if(hit.Charge()>5.0) {
      obj.FillHistogram("uml_summary",6400,0,64000,hit.Charge(),
          40,0,40,hit.GetChannel());
    }
  }

  if(uml->gamma_time>10 && uml->GetPin1().Timestamp()>10) {
    double delta = uml->GetPin1().Timestamp()-uml->gamma_time;
    //std::cout << "delta    " << delta << std::endl;
    //std::cout << "energy   " << uml->gamma_time << std::endl;
   // std::cout << "pin1     " << uml->GetPin1().Timestamp() << std::endl;
   // std::cout << "gt       " << uml->gamma_time << std::endl;
    obj.FillHistogram("ddas_gamma_time",1000,-500,500,delta,
                                        4000,0,10000,uml->gamma_energy);
    obj.FillHistogram("ddas_gamma_with_particle",10000,0,10000,uml->gamma_energy);

  }


  //return;

  uml->CalParameters();

  if(uml->GetTof() <300) return;
  if(uml->GetBeta() < 0.3)  return;
  if(uml->GetBeta() > 0.9999999)  return;

  std::string histname = "";
  std::string dirname  = "";

 // for(unsigned int i=0;i<uml->Size();i++) {
 //   TUMLHit hit = uml->GetUMLHit(i);
 //   if(hit.Charge()>5.0) {
 //     obj.FillHistogram("uml_summary",6400,0,64000,hit.Charge(),
 //         40,0,40,hit.GetChannel());
 //   }
 // }
  dirname = "pid";  

  //if(e1_e2 && e1_e2->IsInside(uml->GetPin1().GetEnergy(),uml->GetPin2().GetEnergy())) {


    obj.FillHistogram(dirname,"pin1_tof",1000,300,500,uml->GetTof(),
        1000,700,1500,uml->GetPin1().GetEnergy());

    obj.FillHistogram(dirname,"pin1_tac",1000,0,300,uml->GetTac1()/1000.,
        750,0,1500,uml->GetPin1().GetEnergy());

    obj.FillHistogram(dirname,"pin2_tof",1000,300,500,uml->GetTof(),
        1000,700,1500,uml->GetPin2().GetEnergy());

    obj.FillHistogram(dirname,"tke_tof",1000,300,500,uml->GetTof(),
        1000,10000,17000,uml->GetTKE());

    obj.FillHistogram(dirname,"pin1_pin2",750,0,1500,uml->GetPin1().GetEnergy(),
        750,0,1500,uml->GetPin2().GetEnergy());

    obj.FillHistogram(dirname,"pin1_tke",750,0,1500,uml->GetPin1().GetEnergy(),
        1000,7000,17000,uml->GetTKE());

    obj.FillHistogram(dirname,"Z_AoQ",1000,2,3,uml->GetAoQ(),
        1000,0,100,uml->GetZ());

    obj.FillHistogram(dirname,"X_AoQ",1000,2,3,uml->GetAoQ(),
        20,-50,50,uml->GetXPosition());

    obj.FillHistogram(dirname,"Z_ZmQ",1000,-2,8,uml->ZmQ(),
        1000,0,100,uml->GetZ());

    obj.FillHistogram(dirname,"Z_Am3Q",1000,-50,0,uml->Am3Q(),
        1000,0,100,uml->GetZ());
  //}
  //std::cout << "tof:  " << uml->GetTof() << std::endl;
  dirname = "uml";
  obj.FillHistogram(dirname,"SSSD",6400,0,6400,uml->GetSssdEnergy());
  obj.FillHistogram(dirname,"PIN1_cal",1500,0,1500,uml->GetPin1().GetEnergy());
  obj.FillHistogram(dirname,"PIN2_cal",1500,0,1500,uml->GetPin2().GetEnergy());
  obj.FillHistogram(dirname,"IMPLANT_cal",17000,0,17000,uml->GetImplant().GetEnergy());




  obj.FillHistogram(dirname,"PIN1",16000,0,64000,uml->GetPin1().Charge());
  obj.FillHistogram(dirname,"PIN2",16000,0,64000,uml->GetPin2().Charge());
  obj.FillHistogram(dirname,"IMPLANT",16000,0,64000,uml->GetImplant().Charge());

  obj.FillHistogram(dirname,"TKE",1700,0,17000,uml->GetTKE());


  //obj.FillHistogram(dirname,"tac1",2000,0,10000,uml->GetTof(),

  obj.FillHistogram(dirname,"AoQ" ,1000,0,3,uml->GetAoQ());
  obj.FillHistogram(dirname,"Z"   ,1000,0,100,uml->GetZ());
  obj.FillHistogram(dirname,"Q"   ,1000,0,100,uml->GetQ());
  obj.FillHistogram(dirname,"A"   ,2500,0,250,uml->A());
  obj.FillHistogram(dirname,"Brho",1000,0,5,uml->GetBrho());
  obj.FillHistogram(dirname,"Beta",1000,0,1,uml->GetBeta());

  obj.FillHistogram(dirname,"Xposition",100,-50,50,uml->GetXPosition());

 // for(int i=0;i<uml->fSssd.size();i++) {
 //   TUMLHit hit = uml->fSssd.at(i);
 //   if(hit.Charge()>10) 
 //     obj.FillHistogram("strips",20,0,20,hit.GetChannel()-15); // 1 -16?
 //   obj.FillHistogram("strips_aligned",6400,0,64000,hit.GetEnergy(),
 //       20,0,20,hit.GetChannel()-15); // 1 -16?
 // }

  dirname = "diag";

  obj.FillHistogram(dirname,"DT_Pin1_Xfp1",1000,-500,500,
      uml->GetPin1().Timestamp() - uml->GetXfp1().Timestamp());
  obj.FillHistogram(dirname,"DT_Pin1_Xfp2",1000,-500,500,
      uml->GetPin1().Timestamp() - uml->GetXfp2().Timestamp());


  uml_map.insert(std::make_pair(uml->Timestamp(),*uml));

}


void HandleGretina(TRuntimeObjects& obj) {
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank88  *bank88 = obj.GetDetector<TBank88>();
  double xfp_time = -1;
  if(bank88) {
    for(unsigned int i=0;i<bank88->Size();i++) {
      TMode3Hit hit = bank88->GetMode3Hit(i);
      if(hit.GetChannel()==0){
        xfp_time = hit.Timestamp();
//        bank88_map.insert(std::make_pair(xfp_time,hit));
      }
    }

  }


  double sum =0.0;
  std::string dirname  = "gretina";
  for(unsigned int i=0;i<gretina->Size();i++) {
    TGretinaHit hit = gretina->GetGretinaHit(i);

    // create map of gretina hits
    if(hit.GetPad() == 0) 
      gretina_map.insert(std::make_pair(hit.Timestamp(),hit));

    obj.FillHistogram(dirname,"summary",2000,0,4000,hit.GetCoreEnergy(),
        200,0,200,hit.GetCrystalId());
    obj.FillHistogram(dirname,"singles",8000,0,8000,hit.GetCoreEnergy());
    if(bank88 && bank88->Timestamp()>0) {
      obj.FillHistogram(dirname,"singles_w88",8000,0,8000,hit.GetCoreEnergy());
      obj.FillHistogram(dirname,"gretinahits_bank88_energy",1000,-500,500,bank88->Timestamp()-hit.GetTime(),
                                                     2000,0,2000,hit.GetCoreEnergy());
      obj.FillHistogram(dirname,"gretinahits_bank88_time",3600,0,3600,hit.Timestamp()/1e8,
                                                          1000,-500,500,bank88->Timestamp()-hit.GetTime());

      if(xfp_time>0) {
        obj.FillHistogram(dirname,"gretinahits_bank88_xfp_energy",1000,-500,500,xfp_time-hit.GetTime(),
                                                                  2000,0,2000,hit.GetCoreEnergy());
        obj.FillHistogram(dirname,"gretinahits_bank88_xfp",3600,0,3600,hit.Timestamp()/1e8,
                                                            1000,-500,500,xfp_time-hit.GetTime());
      }
    } else {
      obj.FillHistogram(dirname,"singles_no88",8000,0,8000,hit.GetCoreEnergy());
    }
    sum += hit.GetCoreEnergy();
  } 
  obj.FillHistogram(dirname,"gretina_sum",8000,0,4000,sum);

}


void HandleDDAS_GRETINA(TRuntimeObjects& obj) {

  //uml->CalParameters();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TUML     *uml     = obj.GetDetector<TUML>();

  if(uml->GetTof() <300) return;
  if(uml->GetBeta() < 0.3)  return;
  if(uml->GetBeta() > 0.9999999)  return;


  double sum =0.0;
  std::string dirname  = "gretinai_particle";
  for(unsigned int i=0;i<gretina->Size();i++) {
    TGretinaHit hit = gretina->GetGretinaHit(i);
    obj.FillHistogram(dirname,"summary",2000,0,4000,hit.GetCoreEnergy(),
        200,0,200,hit.GetCrystalId());
    sum += hit.GetCoreEnergy();
  } 
  obj.FillHistogram(dirname,"gretina_sum",8000,0,4000,sum);




  TList *gates = &(obj.GetGates());
  TIter iter(gates);
  while(TObject *key = iter.Next()) {
    GCutG *gate = (GCutG*)key;
    if(gate->IsInside(uml->GetTof(),uml->GetTKE())) {
      for(unsigned int i=0;i<gretina->Size();i++) {
        TGretinaHit hit = gretina->GetGretinaHit(i);
        obj.FillHistogram(dirname,Form("gamma_%s",gate->GetName()),8000,0,8000,hit.GetCoreEnergy());

        //double deltatime = uml->GetImplant().Timestamp() - hit.Timestamp();
        double deltatime = uml->Timestamp() - hit.Timestamp();
        //std::cout << "delta  " << deltatime << std::endl;
        obj.FillHistogram(dirname,Form("dt_%s",gate->GetName()),500,-3000,3000,deltatime,
                                                                4000,0,8000,hit.GetCoreEnergy());

      }
    }
  }

}



void SearchIosmer(TRuntimeObjects &obj){
  long twin = 3000; // ticks,  30us
  if(gretina_map.size()==0 || bank88_map.size() == 0) return;

  std::vector<std::multimap<long, TGretinaHit>::iterator> store_ghit;
  auto itr = gretina_map.rbegin();
  while ( gretina_map.size()> 0 && bank88_map.size()>0 && (itr->first - bank88_map.begin()->first > twin)) {
    // first, remove gretina events older than bank88 ts - twin
    while(bank88_map.begin()->first - gretina_map.begin()->first > twin){
      gretina_map.erase(gretina_map.begin()) ;
    }

    long fts = bank88_map.begin()->first;

    for( auto it = gretina_map.begin(); it!=gretina_map.end(); it++){
      if(it->first>fts+twin) break;
      TGretinaHit tmp_hit = it->second;
      obj.FillHistogram("gamma_time_bank88",300,-3000,3000,it->first-fts,4000,0,4000,tmp_hit.GetCoreEnergy());
      store_ghit.push_back(it);
    }

    bool *fill = new bool[store_ghit.size()];
    memset(fill,0,sizeof(bool)*(store_ghit.size()));

    if(store_ghit.size()>1)
      for(size_t l = 0; l<store_ghit.size()-1;l++){
        if(store_ghit.at(l+1)->first-store_ghit.at(l)->first < 50){
          if(fill[l] == false){
            TGretinaHit tmp_hit = store_ghit.at(l)->second;
            obj.FillHistogram("gamma_time_bank88_m2",3000,-3000,3000,tmp_hit.Timestamp()-fts,4000,0,4000,tmp_hit.GetCoreEnergy());
            fill[l] = true;
          }
          if(fill[l+1] == false){
            TGretinaHit tmp_hit = store_ghit.at(l+1)->second;
            obj.FillHistogram("gamma_time_bank88_m2",3000,-3000,3000,tmp_hit.Timestamp()-fts,4000,0,4000,tmp_hit.GetCoreEnergy());
            fill[l+1] = true;
          }
        }
      }

    delete [] fill;
    store_ghit.clear();
    bank88_map.erase(bank88_map.begin()); // done with the oldest hit in bank 88
  }

}


void SearchIosmer2(TRuntimeObjects &obj){
  long twin = 3000; // ticks,  30us
  if(gretina_map.size()==0 || uml_map.size() == 0) return;

  std::vector<std::multimap<long, TGretinaHit>::iterator> store_ghit;
  auto itr = gretina_map.rbegin();
  while ( gretina_map.size()> 0 && uml_map.size()>0 && (itr->first - uml_map.begin()->first > twin)) {
    // first, remove gretina events older than bank88 ts - twin
    while(uml_map.begin()->first - gretina_map.begin()->first > twin){
      gretina_map.erase(gretina_map.begin()) ;
    }

    long fts = uml_map.begin()->first;

    for( auto it = gretina_map.begin(); it!=gretina_map.end(); it++){
      if(it->first>fts+twin) break;
      TGretinaHit tmp_hit = it->second;
      obj.FillHistogram("gamma_time_uml",300,-twin,twin,it->first-fts,2000,0,4000,tmp_hit.GetCoreEnergy());
      store_ghit.push_back(it);
      

      // apply pid gates
       TList *gates = &(obj.GetGates());
       TIter iter(gates);
       TUML uml = uml_map.begin()->second;
       while(TObject *key = iter.Next()) {
         GCutG *gate = (GCutG*)key;
         if(gate->IsInside(uml.GetTof(),uml.GetTKE())) {
             obj.FillHistogram(Form("gamma_energy_uml_%s",gate->GetName()),8000,0,8000,tmp_hit.GetCoreEnergy());
     
             double deltatime = tmp_hit.Timestamp()-uml.Timestamp();
             obj.FillHistogram(Form("gamma_time_uml_%s",gate->GetName()),500,-twin,twin,deltatime,
                                                                     4000,0,8000,tmp_hit.GetCoreEnergy());
     
         }
       }


    }

    bool *fill = new bool[store_ghit.size()];
    memset(fill,0,sizeof(bool)*(store_ghit.size()));

    if(store_ghit.size()>1)
      for(size_t l = 0; l<store_ghit.size()-1;l++){
        if(store_ghit.at(l+1)->first-store_ghit.at(l)->first < 50){
          if(fill[l] == false){
            TGretinaHit tmp_hit = store_ghit.at(l)->second;
            obj.FillHistogram("gamma_time_uml_m2",3000,-twin,twin,tmp_hit.Timestamp()-fts,2000,0,4000,tmp_hit.GetCoreEnergy());
            fill[l] = true;
          }
          if(fill[l+1] == false){
            TGretinaHit tmp_hit = store_ghit.at(l+1)->second;
            obj.FillHistogram("gamma_time_uml_m2",3000,-twin,twin,tmp_hit.Timestamp()-fts,2000,0,4000,tmp_hit.GetCoreEnergy());
            fill[l+1] = true;
          }
        }
      }

    delete [] fill;
    store_ghit.clear();
    uml_map.erase(uml_map.begin()); // done with the oldest hit in bank 88
  }

}





// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank88  *bank88  = obj.GetDetector<TBank88>();
  TUML     *uml     = obj.GetDetector<TUML>();

// if(!e1_e2) {
//   TList *gates = &(obj.GetGates());
//   TIter iter(gates);
//   while(TObject *key = iter.Next()) {
//     GCutG *gate = (GCutG*)key;
//     if(!strcmp(gate->GetName(),"e1_e2")) {
//       printf("loaded e1_e2...\n");
//       e1_e2 = gate;
//     }
//   }
// }

  TList    *list    = &(obj.GetObjects());
  int numobj        = list->GetSize();

  std::string histname = "";
  std::string dirname  = "";

  if(uml) {
    HandleUML(obj);
  }

  if(gretina) {
    HandleGretina(obj);
  }

  if(bank88 && uml) {
    obj.FillHistogram("bank88_ddas",3600,0,3600,(bank88->Timestamp())/1e8,
        500,-500,500,bank88->Timestamp()-  uml->Timestamp()  );
  }

  if(uml && gretina) {
    obj.FillHistogram("gretina_ddas",3600,0,3600,(gretina->Timestamp())/1e8,
        500,-500,500,gretina->Timestamp()-  uml->Timestamp());

    obj.FillHistogram("gretina_ddaspin1",17000,0,17000,uml->GetPin1().GetEnergy(),
        500,-500,500,gretina->Timestamp()-  uml->Timestamp());
    obj.FillHistogram("gretina_ddaspin2",17000,0,17000,uml->GetPin2().GetEnergy(),
        500,-500,500,gretina->Timestamp()-  uml->Timestamp());

    for(unsigned int i=0;i<uml->Size();i++) {

      obj.FillHistogram("gretina_ddasihit",40,0,40,uml->GetUMLHit(i).GetChannel(),
          500,-500,500,gretina->Timestamp()-  uml->Timestamp());

    }

    HandleDDAS_GRETINA(obj);
  }


  SearchIosmer2(obj);



  if(numobj!=list->GetSize())
    list->Sort();

}

