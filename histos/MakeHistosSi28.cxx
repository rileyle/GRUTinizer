
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
#include "TBank29.h"
#include "TS800.h"
#include "TS800Sim.h"
#include "TGretSim.h"
#include "GValue.h"


#include "TChannel.h"
#include "GValue.h"

// Relativistic gamma with S800 dta correction
// (following Dirk's SpecTcl GretaCalculator class)
Double_t Gamma(Double_t beta, Double_t dta){

  Double_t gamma =  1/sqrt(1. - beta*beta);

  // DTA correction
  Double_t b = beta;
  Double_t dp_p = gamma/(1. + gamma)*dta;
  b = b*(1. + dp_p /(gamma*gamma));

  return 1./sqrt(1. - b*b);
}

// cos(theta) for Doppler reconstruction with "Phi" correction
// (following Dirk's SpecTcl GretaCalculator class)
Double_t CosTheta(Double_t x, Double_t y, Double_t z,
		  Double_t ata, Double_t bta){

  TVector3 coor(x, y, z);

  // "Phi" correction (correct for incoming ata, bta)
  Double_t xsin = sin(ata); // Offsets are handled by TS800 with GValues
  Double_t ysin = sin(bta); //         ATA_SHIFT and BTA_SHIFT

  Double_t s800_theta;
  Double_t s800_phi;
  if (xsin > 0 && ysin > 0) {
    //s800_phi = -atan(ysin/xsin) ;
    s800_phi = 2*M_PI-atan(ysin/xsin) ;
  } else if (xsin < 0 && ysin > 0) {
    //s800_phi = -M_PI + atan(ysin/fabs(xsin)) ;
    s800_phi = M_PI + atan(ysin/fabs(xsin)) ;
  } else if (xsin < 0 && ysin < 0) {
    //s800_phi = -M_PI - atan(fabs(ysin)/fabs(xsin));
    s800_phi = M_PI - atan(fabs(ysin)/fabs(xsin));
  } else if (xsin > 0 && ysin < 0) {
    //s800_phi = -2*M_PI + atan(fabs(ysin)/xsin);
    s800_phi = atan(fabs(ysin)/xsin);
  } else {
    s800_phi = 0.0 ;
  }
  s800_theta = asin(sqrt(xsin*xsin + ysin*ysin));

  return (sin(coor.Theta())*sin(s800_theta) * 
	  (sin(coor.Phi())*sin(s800_phi)
	   +cos(coor.Phi())*cos(s800_phi))
	  +cos(coor.Theta())*cos(s800_theta));  

}

// cos(theta) for Doppler reconstruction with "Phi" correction
// (following Dirk's SpecTcl GretaCalculator class)
// using the exact direction from the TGretSimHit 
Double_t CosThetaAng(Double_t theta, Double_t phi,
		     Double_t ata, Double_t bta){

  // "Phi" correction (correct for incoming ata, bta)
  Double_t xsin = sin(ata); // Offsets are handled by TS800 with GValues
  Double_t ysin = sin(bta); //         ATA_SHIFT and BTA_SHIFT

  Double_t s800_theta;
  Double_t s800_phi;
  if (xsin > 0 && ysin > 0) {
    //s800_phi = -atan(ysin/xsin) ;
    s800_phi = 2*M_PI-atan(ysin/xsin) ;
  } else if (xsin < 0 && ysin > 0) {
    //s800_phi = -M_PI + atan(ysin/fabs(xsin)) ;
    s800_phi = M_PI + atan(ysin/fabs(xsin)) ;
  } else if (xsin < 0 && ysin < 0) {
    //s800_phi = -M_PI - atan(fabs(ysin)/fabs(xsin));
    s800_phi = M_PI - atan(fabs(ysin)/fabs(xsin));
  } else if (xsin > 0 && ysin < 0) {
    //s800_phi = -2*M_PI + atan(fabs(ysin)/xsin);
    s800_phi = atan(fabs(ysin)/xsin);
  } else {
    s800_phi = 0.0 ;
  }
  s800_theta = asin(sqrt(xsin*xsin + ysin*ysin));

  return (sin(theta)*sin(s800_theta) * 
	  (sin(phi)*sin(s800_phi)
	   +cos(phi)*cos(s800_phi))
	  +cos(theta)*cos(s800_theta));  

}

#define Q1 15
#define Q2 7
#define Q3 8
#define Q4 16
#define Q5 9
#define Q6 14
#define Q7 17
#define Q8 6
#define Q9 19

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
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  TBank29  *bank29  = obj.GetDetector<TBank29>();
  TS800 *s800       = obj.GetDetector<TS800>();
  TS800Sim *s800Sim = obj.GetDetector<TS800Sim>();
  TGretSim *gretSim = obj.GetDetector<TGretSim>();

  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();
  
  if(!s800Sim)
    return;
  
  Double_t dta = s800Sim->GetS800SimHit(0).GetDTA();

  // Rough dta acceptance cut
  if(dta < -0.06 || dta > 0.06)
    return;
  
  obj.FillHistogram("s800","dta",
		    4096, -6, 6,
		    dta*100.);

  Double_t angRes = GValue::Value("S800_ANG_RES");
  if(std::isnan(angRes))
    angRes = 0.;
  Double_t ata = s800Sim->GetS800SimHit(0).GetATA()/1000.
    + gRandom->Gaus(0., angRes);
  Double_t bta = s800Sim->GetS800SimHit(0).GetBTA()/1000.
    + gRandom->Gaus(0., angRes);

  Double_t xsin =  sin(ata);
  Double_t ysin = -sin(bta);
  Double_t scatter = asin(sqrt(xsin*xsin + ysin*ysin))*1000.;

  obj.FillHistogram("s800","scatter",
		    4096, 0, 300,
		    scatter);
    
  if(gretSim){
    obj.FillHistogram("sim","beta",
		      500, 0, 0.5,
		      gretSim->GetGretinaSimHit(0).GetBeta());
    obj.FillHistogram("sim","z",
		      1000,-5,5.,
		      gretSim->GetGretinaSimHit(0).GetZ());

    obj.FillHistogram("sim","beta_z",
		      1000,-5,5.,
		      gretSim->GetGretinaSimHit(0).GetZ(),
		      300, 0.2, 0.5,
		      gretSim->GetGretinaSimHit(0).GetBeta());

    obj.FillHistogram("sim","x",
		      400,-20., 20.,
		      gretSim->GetGretinaSimHit(0).GetX());

    obj.FillHistogram("sim","y",
		      400,-20., 20.,
		      gretSim->GetGretinaSimHit(0).GetY());

    obj.FillHistogram("sim","xy",
		      400,-20., 20.,
		      gretSim->GetGretinaSimHit(0).GetX(),
		      400,-20., 20.,
		      gretSim->GetGretinaSimHit(0).GetY());
    
  }

  if(!gretina)
    return;

  double xoffset = GValue::Value("GRETINA_X_OFFSET");
  if(std::isnan(xoffset))
    xoffset=0.00;
  double yoffset = GValue::Value("GRETINA_Y_OFFSET");
  if(std::isnan(yoffset))
    yoffset=0.00;
  double zoffset = GValue::Value("GRETINA_Z_OFFSET");
  if(std::isnan(zoffset))
    zoffset=0.00;
  TVector3 targetOffset(xoffset,yoffset,zoffset);

  // GRETINA position resolution
  Double_t posRes = GValue::Value("GRETINA_POS_RES");
  if(std::isnan(posRes))
    posRes=0.00;

  // GRETINA energy resolution (intrinsic)
  Double_t eRes = GValue::Value("GRETINA_ENERGY_RES");
  if(std::isnan(eRes))
    eRes=0.0;

  Double_t sx = gRandom->Gaus(0., posRes);
  Double_t sy = gRandom->Gaus(0., posRes);
  Double_t sz = gRandom->Gaus(0., posRes);
  TVector3 positionSmear(sx, sy, sz);
  
  const Int_t    nBetas = 1;
  Double_t betas[nBetas] = {0.378}; // Si28 mid-target
  
  Int_t    energyNChannels = 8192;
  Double_t energyLlim = 0.;
  Double_t energyUlim = 8192.;
  
  Double_t calorimeterEnergy = 0.;
  std::vector<TGretinaHit> hits;
  
  for(int x=0; x<gretina->Size(); x++){

    TGretinaHit hit = gretina->GetGretinaHit(x);

    obj.FillHistogram("position", "theta",
		      180, 0., 180.,
		      hit.GetTheta()*TMath::RadToDeg());

    Double_t xtal_theta = hit.GetCrystalPosition().Theta()*TMath::RadToDeg();
    obj.FillHistogram("position", "xtal_theta",
		      180, 0., 180.,
		      xtal_theta);
    Int_t ring = -1;
    if(xtal_theta < 52.){
      obj.FillHistogram("position", "ring0_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 0;
    } else if(xtal_theta < 60.){
      obj.FillHistogram("position", "ring1_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 1;
    } else if(xtal_theta < 66.){
      obj.FillHistogram("position", "ring2_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 2;
    } else if(xtal_theta < 70.){
      obj.FillHistogram("position", "ring3_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 3;
    } else if(xtal_theta < 81.){
      obj.FillHistogram("position", "ring4_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 4;
    } else if(xtal_theta < 90.){
      obj.FillHistogram("position", "ring5_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 5;
    } else if(xtal_theta < 98.){
      obj.FillHistogram("position", "ring6_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 6;      
    } else {
      obj.FillHistogram("position", "ring7_theta",
			180, 0., 180.,
			xtal_theta);
      ring = 7;
    }
    
    obj.FillHistogram("position", "theta_vs_phi",
		      360, 0., 360.,
		      hit.GetPhi()*TMath::RadToDeg(),
		      180, 0., 180.,
		      hit.GetTheta()*TMath::RadToDeg());

    if(hit.GetHoleNumber() < 10){
      obj.FillHistogram("position", "theta_vs_phi_fw",
			360, 0., 360.,
			hit.GetPhi()*TMath::RadToDeg(),
			180, 0., 180.,
			hit.GetTheta()*TMath::RadToDeg());
      if(hit.GetCrystalNumber() == 2){
	obj.FillHistogram("position", "theta_vs_phi_xfw",
			  360, 0., 360.,
			  hit.GetPhi()*TMath::RadToDeg(),
			  180, 0., 180.,
			  hit.GetTheta()*TMath::RadToDeg());
      }
    } else {
      obj.FillHistogram("position", "theta_vs_phi_bw",
			360, 0., 360.,
			hit.GetPhi()*TMath::RadToDeg(),
			180, 0., 180.,
			hit.GetTheta()*TMath::RadToDeg());
      if(hit.GetHoleNumber() > 13 && hit.GetCrystalNumber() == 1){
	obj.FillHistogram("position", "theta_vs_phi_xbw",
			  360, 0., 360.,
			  hit.GetPhi()*TMath::RadToDeg(),
			  180, 0., 180.,
			  hit.GetTheta()*TMath::RadToDeg());
      }
    }

    Double_t coreEnergy = hit.GetCoreEnergy()*gRandom->Gaus(1., eRes);

    // Full energy event?
    Bool_t fullEnergy = false;
    Double_t de = gretSim->GetGretinaSimHit(0).GetEn()
      - hit.GetCoreEnergy();
    if( fabs(de) < 0.1 )
      fullEnergy = true;

    obj.FillHistogram("energy",  "de",
		      200, -10., 10.,
		      de);
    
    // Position of main interaction (largest energy deposit)
    TVector3 mainIntPos = hit.GetInteractionPosition(hit.GetFirstIntPoint())
                        - targetOffset + positionSmear;

    Double_t xInt = mainIntPos.X();
    Double_t yInt = mainIntPos.Y();
    Double_t zInt = mainIntPos.Z();

    // Exact direction of gamma emission
    
    Double_t thSim = gretSim->GetGretinaSimHit(0).GetTheta();
    Double_t phSim = gretSim->GetGretinaSimHit(0).GetPhi();//+TMath::PiOver2();

    // "Measured" theta of the first interaction
    Int_t iFirstHit = -1;
    TVector3 firstIntPos;
    Int_t iZerothHit = -1;
    TVector3 zerothIntPos;

    // std::cout << "New hit" << std::endl;
    // std::cout << "  IPs:" << std::endl;
    for(Int_t ip = 0; ip < hit.NumberOfInteractions(); ip++){

      TVector3 ipPos = hit.GetInteractionPosition(ip);
      TVector3 emPos = TVector3(gretSim->GetGretinaSimHit(0).GetX(),
				gretSim->GetGretinaSimHit(0).GetY(),
				gretSim->GetGretinaSimHit(0).GetZ());
      TVector3 correctedPos = ipPos - emPos;

      // std::cout << "   " << ip << "  " << ipPos.X() << ", "
      // 		<< ipPos.Y() << ", "
      // 		<< ipPos.Z() << ", e = "
      // 		<< hit.GetSegmentEng(ip) << std::endl;

      if( fabs(correctedPos.Theta() - thSim) < 0.008 )
	iFirstHit = ip;

      //      if( fabs(hit.GetZerothInteractionEnergy()
      //	       - hit.GetSegmentEng(ip)) < 0.001 ) {
      if( hit.GetZerothInteractionPosition() == ipPos ){
        if(iZerothHit > 0)
	  std::cout << "**** Holly crap! identical positions!!!! ****"
		    << std::endl;

	// else
	//   iZerothHit = ip; // Take the first one

	iZerothHit = ip; // Take the last one
      }
    }

    // std::cout << "  main: 0  " << xInt << ", " << yInt << ", " << zInt
    // 	      << std::endl;

    // TVector3 tmpPos = hit.GetInteractionPosition(iZerothHit);
    // std::cout << "   0th: " << iZerothHit << "  "
    // 	      << tmpPos.X() << ", "
    // 	      << tmpPos.Y() << ", "
    // 	      << tmpPos.Z()
    // 	      << std::endl;
    // tmpPos = hit.GetInteractionPosition(iFirstHit);
    // std::cout << "   1st: " << iFirstHit << "  "
    // 	      << tmpPos.X() << ", "
    // 	      << tmpPos.Y() << ", "
    // 	      << tmpPos.Z()
    // 	      << std::endl;

    if(iFirstHit == -1){
      if(fullEnergy){
	std::cout << "*** Failed to find a first hit" << std::endl;
      }
    } else {
      firstIntPos = hit.GetInteractionPosition(iFirstHit)
                  - targetOffset + positionSmear;
    }
    
    if(iZerothHit == -1){
      if(fullEnergy){
	std::cout << "*** Failed to find a zeroth hit" << std::endl;
      }
    } else {
      zerothIntPos = hit.GetInteractionPosition(iZerothHit)
                   - targetOffset + positionSmear;
    }
    
    Double_t xInt1st = firstIntPos.X();
    Double_t yInt1st = firstIntPos.Y();
    Double_t zInt1st = firstIntPos.Z();

    Double_t xInt0th = zerothIntPos.X();
    Double_t yInt0th = zerothIntPos.Y();
    Double_t zInt0th = zerothIntPos.Z();
    
    //    thSim += gRandom->Gaus(0., posRes/185.);  // Simulate GT resolution
    //    thSim += gRandom->Gaus(0., angRes/1000.); // Simulate S800 resolution

    //    phSim += gRandom->Gaus(0., posRes/185.);  // Simulate GT resolution
    //    phSim += gRandom->Gaus(0., angRes/1000.); // Simulate S800 resolution

    // Addback preprocessing
    if(coreEnergy > energyLlim &&
       coreEnergy < energyUlim){

      calorimeterEnergy      += coreEnergy;

      hits.push_back(hit);

    }
    
    //                directory, histogram
    obj.FillHistogram("energy",  "overview",
		      energyNChannels, energyLlim, energyUlim,
		      coreEnergy,
		      100, 0, 100, hit.GetCrystalId());

    obj.FillHistogram("energy",  "energy",
		      energyNChannels, energyLlim, energyUlim,
		      coreEnergy);

    obj.FillHistogram("position",
		      "thMain",
		      1800, 0., 180.,
		      mainIntPos.Theta()*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "phiMain",
		      3600, -180., 180.,
		      mainIntPos.Phi()*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "th0th",
		      1800, 0., 180.,
		      zerothIntPos.Theta()*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "phi0th",
		      3600, -180., 180.,
		      zerothIntPos.Phi()*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "thSim",
		      1800, 0., 180.,
		      thSim*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "phSim",
		      3600, -180., 180.,
		      phSim*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dthetaMain",
		      4000, -20., 20.,
		      (mainIntPos.Theta() - thSim)*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dphiMain",
		      4000, -20., 20.,
		      (mainIntPos.Phi() - phSim)*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dtheta0th",
		      4000, -20., 20.,
		      (zerothIntPos.Theta() - thSim)*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dphi0th",
		      4000, -20., 20.,
		      (zerothIntPos.Phi() - phSim)*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dtheta0th_overview",
		      4000, -20., 20.,
		      (zerothIntPos.Theta() - thSim)*TMath::RadToDeg(),
		      100, 0, 100, hit.GetCrystalId());
    obj.FillHistogram("position",
		      "dphi0th_overview",
		      4000, -20., 20.,
		      (zerothIntPos.Phi() - phSim)*TMath::RadToDeg(),
		      100, 0, 100, hit.GetCrystalId());
    obj.FillHistogram("position",
		      "dtheta1st",
		      4000, -20., 20.,
		      (firstIntPos.Theta() - thSim)*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dthetaMain0th",
		      4000, -20., 20.,
		      (mainIntPos.Theta()
		       - zerothIntPos.Theta())*TMath::RadToDeg());
    obj.FillHistogram("position",
		      "dthetaMain1st",
		      4000, -20., 20.,
		      (mainIntPos.Theta()
		       - firstIntPos.Theta())*TMath::RadToDeg());
    if(fullEnergy){
      obj.FillHistogram("position",
			"thMain_full",
			1800, 0., 180.,
			mainIntPos.Theta()*TMath::RadToDeg());
      obj.FillHistogram("position",
			"phiMain_full",
			3600, -180., 180.,
			mainIntPos.Phi()*TMath::RadToDeg());
      obj.FillHistogram("position",
			"th0th_full",
			1800, 0., 180.,
			zerothIntPos.Theta()*TMath::RadToDeg());
      obj.FillHistogram("position",
			"phi0th_full",
			3600, -180., 180.,
			zerothIntPos.Phi()*TMath::RadToDeg());
      obj.FillHistogram("position",
			"thSim_full",
			1800, 0., 180.,
			thSim*TMath::RadToDeg());
      obj.FillHistogram("position",
			"phSim_full",
			3600, -180., 180.,
			phSim*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dthetaMain_full",
			4000, -20., 20.,
			(mainIntPos.Theta() - thSim)*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dphiMain_full",
			4000, -20., 20.,
			(mainIntPos.Phi() - phSim)*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dtheta0th_full",
			4000, -20., 20.,
			(zerothIntPos.Theta() - thSim)*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dphi0th_full",
			4000, -20., 20.,
			(zerothIntPos.Phi() - phSim)*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dtheta1st_full",
			4000, -20., 20.,
			(firstIntPos.Theta() - thSim)*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dthetaMain0th_full",
			4000, -20., 20.,
			(mainIntPos.Theta()
			 - zerothIntPos.Theta())*TMath::RadToDeg());
      obj.FillHistogram("position",
			"dthetaMain1st_full",
			4000, -20., 20.,
			(mainIntPos.Theta()
			 - firstIntPos.Theta())*TMath::RadToDeg());
    }
    
    for(int i=0; i<nBetas; i++){

      // Various Doppler corrections
      Double_t gamma0 = 1./sqrt(1. - betas[i]*betas[i]);
      Double_t costheta0 = cos(mainIntPos.Theta());
      
      Double_t gamma  = Gamma(betas[i], dta);
      Double_t costheta  = CosTheta(xInt, yInt, zInt, ata, bta);

      // Double_t costheta1st = CosTheta(xInt1st, yInt1st, zInt1st, ata, bta);

      Double_t costheta0th = CosTheta(xInt0th, yInt0th, zInt0th, ata, bta);

      Double_t costhetasim = CosThetaAng(thSim, phSim, ata, bta);

      Double_t doppler0   = gamma0*(1. - betas[i]* costheta0);
      Double_t dopplerPhi = gamma0*(1. - betas[i]* costheta);
      Double_t dopplerDta = gamma*(1. - betas[i]* costheta0);
      Double_t doppler    = gamma*(1. - betas[i]* costheta);

      // Double_t doppler1st  = gamma*(1. - betas[i]* costheta1st);

      Double_t doppler0th  = gamma*(1. - betas[i]* costheta0th);

      Double_t dopplersim  = gamma*(1. - betas[i]* costhetasim);
      
      Double_t dopEnergy0   = doppler0*coreEnergy;
      Double_t dopEnergyPhi = dopplerPhi*coreEnergy;
      Double_t dopEnergyDta = dopplerDta*coreEnergy;
      Double_t dopEnergy    = doppler*coreEnergy;

      // Double_t dopEnergy1st  = doppler1st*coreEnergy;

      Double_t dopEnergy0th  = doppler0th*coreEnergy;

      Double_t dopEnergySim  = dopplersim*coreEnergy;

      obj.FillHistogram("energy",
			Form("dop0_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy0);
      obj.FillHistogram("energy",
			Form("dopPhi_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergyPhi);
      obj.FillHistogram("energy",
			Form("dopDta_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergyDta);
      obj.FillHistogram("energy",
			Form("dop_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy);
      obj.FillHistogram("energy",
			Form("dop0th_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy0th);
      // obj.FillHistogram("energy",
      // 			Form("dop1st_%.0f", betas[i]*10000),
      // 			energyNChannels, energyLlim, energyUlim,
      // 			dopEnergy1st);
      obj.FillHistogram("energy",
			Form("dopSim_%.0f", betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergySim);

      // obj.FillHistogram("energy",
      // 			Form("dopPhi_vs_thetaMain_%.0f",
      // 			     betas[i]*10000),
      // 			180, 0., 180.,
      // 			mainIntPos.Theta()*TMath::RadToDeg(),
      // 			energyNChannels, energyLlim, energyUlim,
      // 			dopEnergy);
      // obj.FillHistogram("energy",
      // 			Form("dopPhi_vs_dthetaMain_%.0f",
      // 			     betas[i]*10000),
      // 			400, -20., 20.,
      // 			(mainIntPos.Theta() - thSim)*TMath::RadToDeg(),
      // 			energyNChannels, energyLlim, energyUlim,
      // 			dopEnergy);

      if(fullEnergy){
	obj.FillHistogram("energy",
			  Form("dop0_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0);
	obj.FillHistogram("energy",
			  Form("dopPhi_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyPhi);
	obj.FillHistogram("energy",
			  Form("dopDta_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyDta);
	obj.FillHistogram("energy",
			  Form("dop_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy);
	obj.FillHistogram("energy",
			  Form("dop0th_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0th);
	// obj.FillHistogram("energy",
	// 		  Form("dop1st_full_%.0f", betas[i]*10000),
	// 		  energyNChannels, energyLlim, energyUlim,
	// 		  dopEnergy1st);
	obj.FillHistogram("energy",
			  Form("dopSim_full_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergySim);

	obj.FillHistogram("energy",
			  Form("dopPhi_vs_theta_full_%.0f",
			       betas[i]*10000),
			  180, 0., 180.,
			  mainIntPos.Theta()*TMath::RadToDeg(),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy);
	obj.FillHistogram("energy",
			  Form("dopPhi_vs_dtheta_full_%.0f",
			       betas[i]*10000),
			  400, -20., 20.,
			  (mainIntPos.Theta() - thSim)*TMath::RadToDeg(),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyPhi);
	obj.FillHistogram("energy",
			  Form("dop0_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0);
	obj.FillHistogram("energy",
			  Form("dopPhi_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyPhi);
	obj.FillHistogram("energy",
			  Form("dopDta_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyDta);
	obj.FillHistogram("energy",
			  Form("dop_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy);
	obj.FillHistogram("energy",
			  Form("dop0th_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0th);
	// obj.FillHistogram("energy",
	// 		  Form("dop1st_ring%d_full_%.0f", ring, betas[i]*10000),
	// 		  energyNChannels, energyLlim, energyUlim,
	// 		  dopEnergy1st);
	obj.FillHistogram("energy",
			  Form("dopSim_ring%d_full_%.0f", ring, betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergySim);
      }
      
      obj.FillHistogram("energy",
			Form("dop0_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy0);
      obj.FillHistogram("energy",
			Form("dopPhi_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergyPhi);
      obj.FillHistogram("energy",
			Form("dopDta_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergyDta);
      obj.FillHistogram("energy",
			Form("dop_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy);
      obj.FillHistogram("energy",
			Form("dop0th_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergy0th);
      // obj.FillHistogram("energy",
      // 			Form("dop1st_ring%d_%.0f", ring, betas[i]*10000),
      // 			energyNChannels, energyLlim, energyUlim,
      // 			dopEnergy1st);
      obj.FillHistogram("energy",
			Form("dopSim_ring%d_%.0f", ring, betas[i]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopEnergySim);
      
      if(hit.GetHoleNumber() < 10){
	obj.FillHistogram("energy",
			  Form("dop0_fw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0);
	obj.FillHistogram("energy",
			  Form("dopPhi_fw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyPhi);
	obj.FillHistogram("energy",
			  Form("dopDta_fw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyDta);
	obj.FillHistogram("energy",
			  Form("dop_fw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy);
	if(hit.GetCrystalNumber() == 2){
	  obj.FillHistogram("energy",
			    Form("dop0_xfw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergy0);
	  obj.FillHistogram("energy",
			    Form("dopPhi_xfw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergyPhi);
	  obj.FillHistogram("energy",
			    Form("dopDta_xfw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergyDta);
	  obj.FillHistogram("energy",
			    Form("dop_xfw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergy);
	}
      } else {
	obj.FillHistogram("energy",
			  Form("dop0_bw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy0);
	obj.FillHistogram("energy",
			  Form("dopPhi_bw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyPhi);
	obj.FillHistogram("energy",
			  Form("dopDta_bw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergyDta);
	obj.FillHistogram("energy",
			  Form("dop_bw_%.0f", betas[i]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopEnergy);
	if(hit.GetHoleNumber() > 13 && hit.GetCrystalNumber() == 1){
	  obj.FillHistogram("energy",
			    Form("dop0_xbw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergy0);
	  obj.FillHistogram("energy",
			    Form("dopPhi_xbw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergyPhi);
	  obj.FillHistogram("energy",
			    Form("dopDta_xbw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergyDta);
	  obj.FillHistogram("energy",
			    Form("dop_xbw_%.0f", betas[i]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopEnergy);
	}
      }
    }
        
  }
  
  // Addback
  obj.FillHistogram("addback",  "calorimeter",
		    energyNChannels, energyLlim, energyUlim,
		    calorimeterEnergy);

  while(hits.size() > 0){
    TGretinaHit currentHit = hits.back();
    hits.pop_back();
    
    // Find and add all hits in a cluster of adjacent crystals including
    // the current hit.
    //
    // CAUTION: This clustering includes neighbors of neighbors!
    std::vector<TGretinaHit> cluster;
    cluster.push_back(currentHit);
    int lastClusterSize = 0;
    while(lastClusterSize < cluster.size()){
      for(int i = 0; i < cluster.size(); i++){
	for(int j = 0; j < hits.size(); j++){
	  TVector3 distance = cluster[i].GetCrystalPosition()
	                       - hits[j].GetCrystalPosition();

	  obj.FillHistogram("position",  "crystal_separation",
			    1000, 0., 1000.,
			    distance.Mag());

	  if(distance.Mag() < 80.){ // Neighbors
	    cluster.push_back(hits.back());
	    hits.pop_back();
	  }
	}
      }
      lastClusterSize = cluster.size();
    }
    
    // Calculate the total energy deposited in the cluster,
    // and count the pairs of neighbors.
    Int_t neighbors = 0;
    Double_t addbackEnergy = 0.;
    TVector3 firstHitPos;
    Int_t firstHitHoleNum;
    Int_t firstHitCrystal;
    Double_t firstHitEnergy = 0;
    for(int i = 0; i < cluster.size(); i++){
      addbackEnergy += cluster[i].GetCoreEnergy()*gRandom->Gaus(1., eRes);
      
      // Find the largest IP in the cluster and save its position
      // for Doppler correction.
      if(cluster[i].GetSegmentEng(cluster[i].GetFirstIntPoint())
	 > firstHitEnergy){
	firstHitHoleNum = cluster[i].GetHoleNumber();
	firstHitCrystal = cluster[i].GetCrystalNumber();
	firstHitPos = cluster[i].GetInteractionPosition(cluster[i].GetFirstIntPoint()) - targetOffset;
	firstHitEnergy = cluster[i].GetSegmentEng(cluster[i].GetFirstIntPoint());
      }
      
      for(int j = i+1; j < cluster.size(); j++){
	TVector3 distance =   cluster[i].GetCrystalPosition()
	                    - cluster[j].GetCrystalPosition();
	if(distance.Mag() < 80.) neighbors++;
      }
    }

    // GRETINA position resolution
    firstHitPos += positionSmear;
    
    obj.FillHistogram("position",  "firstHitHoleNumber",
		      30, 0, 29, 
		      firstHitHoleNum);
        
    // Doppler correct the addback energy.
    Double_t dopplerABEnergy0[nBetas]   = {0.};
    Double_t dopplerABEnergyPhi[nBetas] = {0.};
    Double_t dopplerABEnergyDta[nBetas] = {0.};
    Double_t dopplerABEnergy[nBetas]    = {0.};
    for(int b=0; b<nBetas; b++){

      Double_t gamma0 = 1./sqrt(1. - betas[b]*betas[b]);
      Double_t costheta0 = cos(firstHitPos.Theta());
      
      Double_t gamma  = Gamma(betas[b], dta);
      Double_t costheta  = CosTheta(firstHitPos.X(),
				    firstHitPos.Y(),
				    firstHitPos.Z(),
				    ata, bta);

      dopplerABEnergy0[b]   = gamma0*(1. - betas[b]*costheta0)*addbackEnergy;
      dopplerABEnergyPhi[b] = gamma0*(1. - betas[b]*costheta)*addbackEnergy;
      dopplerABEnergyDta[b] = gamma*(1. - betas[b]*costheta0)*addbackEnergy;
      dopplerABEnergy[b]    = gamma*(1. - betas[b]*costheta)*addbackEnergy;
    }
    
    TString addbackType;
    if(neighbors == 0 && cluster.size() == 1)
      addbackType = "n0";
    else if(neighbors == 1 && cluster.size() == 2)
      addbackType = "n1";
    else if(neighbors == 3 && cluster.size() == 3)
      addbackType = "n2";
    else
      addbackType = "ng";

    // Fill addback histograms.

    obj.FillHistogram("addback",  addbackType,
		      energyNChannels, energyLlim, energyUlim,
		      addbackEnergy);

    for(int b=0; b<nBetas; b++){
      obj.FillHistogram("addback",
			Form("dop0_%s_%.0f",
			     addbackType.Data(),
			     betas[b]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopplerABEnergy0[b]);
      obj.FillHistogram("addback",
			Form("dopPhi_%s_%.0f",
			     addbackType.Data(),
			     betas[b]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopplerABEnergyPhi[b]);
      obj.FillHistogram("addback",
			Form("dopDta_%s_%.0f",
			     addbackType.Data(),
			     betas[b]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopplerABEnergyDta[b]);
      obj.FillHistogram("addback",
			Form("dop_%s_%.0f",
			     addbackType.Data(),
			     betas[b]*10000),
			energyNChannels, energyLlim, energyUlim,
			dopplerABEnergy[b]);
      if(firstHitHoleNum < 10){
	obj.FillHistogram("addback",
			  Form("dop0_fw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy0[b]);
	obj.FillHistogram("addback",
			  Form("dopPhi_fw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyPhi[b]);
	obj.FillHistogram("addback",
			  Form("dopDta_fw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyDta[b]);
	obj.FillHistogram("addback",
			  Form("dop_fw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy[b]);
	if(firstHitCrystal == 2){
	  obj.FillHistogram("addback",
			    Form("dop0_xfw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_xfw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_xfw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_xfw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	}
      } else {
	obj.FillHistogram("addback",
			  Form("dop0_bw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy0[b]);
	obj.FillHistogram("addback",
			  Form("dopPhi_bw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyPhi[b]);
	obj.FillHistogram("addback",
			  Form("dopDta_bw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyDta[b]);
	obj.FillHistogram("addback",
			  Form("dop_bw_%s_%.0f",
			       addbackType.Data(),
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy[b]);
	if(firstHitHoleNum > 13 && firstHitCrystal == 1){
	  obj.FillHistogram("addback",
			    Form("dop0_xbw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_xbw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_xbw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_xbw_%s_%.0f",
				 addbackType.Data(),
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	}
      }
    }
    
    if(addbackType == "n0"
       || addbackType == "n1"){
      obj.FillHistogram("addback",  "n0n1",
			energyNChannels, energyLlim, energyUlim,
			addbackEnergy);
      for(int b=0; b<nBetas; b++){
	obj.FillHistogram("addback",
			  Form("dop0_n0n1_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy0[b]);
	obj.FillHistogram("addback",
			  Form("dopPhi_n0n1_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyPhi[b]);
	obj.FillHistogram("addback",
			  Form("dopDta_n0n1_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyDta[b]);
	obj.FillHistogram("addback",
			  Form("dop_n0n1_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy[b]);
	// obj.FillHistogram("addback",
	// 		  Form("dop_vs_theta_n0n1_%.0f",
	// 		       betas[b]*10000),
	// 		  180, 0., 180.,
	// 		  firstHitPos.Theta()*TMath::RadToDeg(),
	// 		  energyNChannels, energyLlim, energyUlim,
	// 		  dopplerABEnergy[b]);
	
	if(firstHitHoleNum < 10){
	  obj.FillHistogram("addback",
			    Form("dop0_fw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_fw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_fw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_fw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	  if(firstHitCrystal == 2){
	    obj.FillHistogram("addback",
			      Form("dop0_xfw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy0[b]);
	    obj.FillHistogram("addback",
			      Form("dopPhi_xfw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyPhi[b]);
	    obj.FillHistogram("addback",
			      Form("dopDta_xfw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyDta[b]);
	    obj.FillHistogram("addback",
			      Form("dop_xfw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy[b]);
	  }
	} else {
	  obj.FillHistogram("addback",
			    Form("dop0_bw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_bw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_bw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_bw_n0n1_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	  if(firstHitHoleNum > 13 && firstHitCrystal == 1){
	    obj.FillHistogram("addback",
			      Form("dop0_xbw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy0[b]);
	    obj.FillHistogram("addback",
			      Form("dopPhi_xbw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyPhi[b]);
	    obj.FillHistogram("addback",
			      Form("dopDta_xbw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyDta[b]);
	    obj.FillHistogram("addback",
			      Form("dop_xbw_n0n1_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy[b]);
	  }
	}
      }
    }
    
    if(addbackType == "n0"
       || addbackType == "n1"
       || addbackType == "n2"){
      obj.FillHistogram("addback",  "n0n1n2",
			energyNChannels, energyLlim, energyUlim,
			addbackEnergy);
      for(int b=0; b<nBetas; b++){
	obj.FillHistogram("addback",
			  Form("dop0_n0n1n2_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy0[b]);
	obj.FillHistogram("addback",
			  Form("dopPhi_n0n1n2_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyPhi[b]);
	obj.FillHistogram("addback",
			  Form("dopDta_n0n1n2_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergyDta[b]);
	obj.FillHistogram("addback",
			  Form("dop_n0n1n2_%.0f",
			       betas[b]*10000),
			  energyNChannels, energyLlim, energyUlim,
			  dopplerABEnergy[b]);
	if(firstHitHoleNum < 10){
	  obj.FillHistogram("addback",
			    Form("dop0_fw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_fw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_fw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_fw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	  if(firstHitCrystal == 2){
	    obj.FillHistogram("addback",
			      Form("dop0_xfw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy0[b]);
	    obj.FillHistogram("addback",
			      Form("dopPhi_xfw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyPhi[b]);
	    obj.FillHistogram("addback",
			      Form("dopDta_xfw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyDta[b]);
	    obj.FillHistogram("addback",
			      Form("dop_xfw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy[b]);
	  }
	} else {
	  obj.FillHistogram("addback",
			    Form("dop0_bw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy0[b]);
	  obj.FillHistogram("addback",
			    Form("dopPhi_bw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyPhi[b]);
	  obj.FillHistogram("addback",
			    Form("dopDta_bw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergyDta[b]);
	  obj.FillHistogram("addback",
			    Form("dop_bw_n0n1n2_%.0f",
				 betas[b]*10000),
			    energyNChannels, energyLlim, energyUlim,
			    dopplerABEnergy[b]);
	  if(firstHitHoleNum > 13 && firstHitCrystal == 1){
	    obj.FillHistogram("addback",
			      Form("dop0_xbw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy0[b]);
	    obj.FillHistogram("addback",
			      Form("dopPhi_xbw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyPhi[b]);
	    obj.FillHistogram("addback",
			      Form("dopDta_xbw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergyDta[b]);
	    obj.FillHistogram("addback",
			      Form("dop_xbw_n0n1n2_%.0f",
				   betas[b]*10000),
			      energyNChannels, energyLlim, energyUlim,
			      dopplerABEnergy[b]);
	  }
	}
      }
    }
    
    obj.FillHistogram("addback",  "clusterSize_vs_neighborPairs",
		      20, 0, 20, neighbors,
		      10, 0, 10, cluster.size());
    
  }
  
  if(numobj!=list->GetSize())
    list->Sort();

}
