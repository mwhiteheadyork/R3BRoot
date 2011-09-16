#include "R3BTra2pDigitizer.h"
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"


// includes for modeling
#include "TGeoManager.h"
#include "TParticle.h"
#include "TVirtualMC.h"
#include "TGeoMatrix.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoShapeAssembly.h"


#include "TVector3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include <string>
#include <iostream>


#include "R3BTraPoint.h"
#include "R3BMCTrack.h"
		
using std::cout;
using std::endl;
		

R3BTra2pDigitizer::R3BTra2pDigitizer() :
  FairTask("R3B Tra2p Digitization scheme ") { 
}


R3BTra2pDigitizer::~R3BTra2pDigitizer() {
}


void R3BTra2pDigitizer::SetParContainers() {

  // Get run and runtime database
  FairRunAna* run = FairRunAna::Instance();
  if ( ! run ) Fatal("SetParContainers", "No analysis run");

  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  if ( ! rtdb ) Fatal("SetParContainers", "No runtime database");

  fTra2pDigiPar = (R3BTra2pDigiPar*)(rtdb->getContainer("R3BTra2pDigiPar"));

  if ( fTra2pDigiPar ) {
      cout << "-I- R3BTra2pDigitizer::SetParContainers() "<< endl;
      cout << "-I- Container R3BTra2pDigiPar  loaded " << endl;
  }

}




InitStatus R3BTra2pDigitizer::Init() {

//  cout<<"Init "<<endl;
  // Get input array 
  FairRootManager* ioman = FairRootManager::Instance();
  if ( ! ioman ) Fatal("Init", "No FairRootManager");
  fTra2pPoints = (TClonesArray*) ioman->GetObject("TraPoint");
  fTra2pMCTrack = (TClonesArray*) ioman->GetObject("MCTrack");
  fMCTrack = (TClonesArray*) ioman->GetObject("MCTrack");
  
   
  // Register output array Tra2pDigi
  fTra2pDigi = new TClonesArray("R3BTra2pDigi",1000);
  ioman->Register("Tra2pDigi", "Digital response in Tra2p", fTra2pDigi, kTRUE);
  
  eventNoTra=0;
  
  // Initialise control histograms


  return kSUCCESS;

}


void R3BTra2pDigitizer::Exec(Option_t* opt) {

   Reset();
   eventNoTra+=1;
//     if(eventNoTra/1000. == (int)eventNoTra/1000.) cout<<"Event #: "<<eventNoTra-1<<endl;
     
     Int_t nentriesTra = fTra2pPoints->GetEntries();
     Int_t nentries = fMCTrack->GetEntries();
     
     Int_t TrackIdTra=0;
     Int_t TrackId=0;
     
     Double_t ss03_se_p1;
     Double_t ss03_spos_p1;
     Double_t ss03_sbw_p1;
     Double_t ss03_sarea_p1;
     Double_t ss03_seta_p1;

     Double_t ss03_ke_p1;
     Double_t ss03_kpos_p1;
     Double_t ss03_kbw_p1;
     Double_t ss03_karea_p1;
     Double_t ss03_keta_p1;

     Double_t ss06_se_p1;
     Double_t ss06_spos_p1;
     Double_t ss06_sbw_p1;
     Double_t ss06_sarea_p1;
     Double_t ss06_seta_p1;

     Double_t ss06_ke_p1;
     Double_t ss06_kpos_p1;
     Double_t ss06_kbw_p1;
     Double_t ss06_karea_p1;
     Double_t ss06_keta_p1;
     
     Double_t ss03_se_p2;
     Double_t ss03_spos_p2;
     Double_t ss03_sbw_p2;
     Double_t ss03_sarea_p2;
     Double_t ss03_seta_p2;

     Double_t ss03_ke_p2;
     Double_t ss03_kpos_p2;
     Double_t ss03_kbw_p2;
     Double_t ss03_karea_p2;
     Double_t ss03_keta_p2;

     Double_t ss06_se_p2;
     Double_t ss06_spos_p2;
     Double_t ss06_sbw_p2;
     Double_t ss06_sarea_p2;
     Double_t ss06_seta_p2;

     Double_t ss06_ke_p2;
     Double_t ss06_kpos_p2;
     Double_t ss06_kbw_p2;
     Double_t ss06_karea_p2;
     Double_t ss06_keta_p2;
     


//******************** SSTs **************************//   


   for (Int_t l=0;l<nentriesTra;l++){
//   cout<<"entries "<<l<<endl;
     
     R3BTraPoint *Tra2p_obj = (R3BTraPoint*) fTra2pPoints->At(l);

     Int_t DetID = Tra2p_obj->GetDetectorID();
     Double_t fX_In = Tra2p_obj->GetXIn();
     Double_t fY_In = Tra2p_obj->GetYIn();
     Double_t fX_Out = Tra2p_obj->GetXOut();
     Double_t fY_Out = Tra2p_obj->GetYOut();
     TrackIdTra = Tra2p_obj->GetTrackID();
     R3BMCTrack *aTrack = (R3BMCTrack*) fTra2pMCTrack->At(TrackIdTra);   
     Int_t PID = aTrack->GetPdgCode();
     Int_t mother = aTrack->GetMotherId();

     Double_t fX = ((fX_In + fX_Out)/2);
     Double_t fY = ((fY_In + fY_Out)/2);     

     
     if(PID==2212 && mother<0){
     
      if (DetID==15 && l==0)
      {
        ss03_se_p1 = 50;
//	ss03_spos_p1 = (-0.038230 + (3.5 + fX)); //detectors couldn't be shifted in simulation, they are shifted here 0.038230
	ss03_spos_p1 = (3.5 + fX);
//	ss03_spos_p1 = (-0.039 + (3.5 + fX));  //Christoph
	ss03_sbw_p1 = 2;
	ss03_sarea_p1 = 0;
	ss03_seta_p1 = 0;
	  
	ss03_ke_p1 = 50;
//	ss03_kpos_p1 = (3.9936-((2 + fY)-0.006402)); //detectors couldn't be shifted in simulation, they are shifted here -0.006402
//	ss03_kpos_p1 = (-0.006402 + (2 + fY));
        ss03_kpos_p1 = (2 + fY);
//	ss03_kpos_p1 = ((-0.01035 + (2 + fY)));  //Christoph (no swop)!!!
	ss03_kbw_p1 = 2;   //kpos for 1 SST swop (with 3.9936-) - requirement for tracker
	ss03_karea_p1 = 0;
	ss03_keta_p1 = 0;
//	cout<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID<<endl;
      }    
      
      if (DetID==15 && l>1)
      {
        ss03_se_p2 = 50;
//	ss03_spos_p2 = (-0.038230 + (3.5 + fX)); //detectors couldn't be shifted in simulation, they are shifted here 0.038230
	ss03_spos_p2 = (3.5 + fX);
//	ss03_spos_p2 = (-0.039 + (3.5 + fX));  //Christoph
	ss03_sbw_p2 = 2;
	ss03_sarea_p2 = 0;
	ss03_seta_p2 = 0;
	  
	ss03_ke_p2 = 50;
//	ss03_kpos_p2 = (3.9936-((2 + fY)-0.006402)); //detectors couldn't be shifted in simulation, they are shifted here -0.006402
//	ss03_kpos_p2 = (-0.006402 + (2 + fY));
	ss03_kpos_p2 = (2 + fY);
//	ss03_kpos_p2 = ((-0.01035 + (2 + fY)));  //Christoph (no swop)!!!
	ss03_kbw_p2 = 2;  //kpos for 1 SST swop (with 3.9936-) - requirement for tracker
	ss03_karea_p2 = 0;
	ss03_keta_p2 = 0;
//	cout<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID<<endl;
      }
            
      if (DetID==16 && l==1)
      {
        ss06_se_p1 = 50;
//	ss06_spos_p1 = (0.038495 + (3.5 + fX));//detectors couldn't be shifted in simulation, they are shifted here -0.038495
	ss06_spos_p1 = (3.5 + fX);
//	ss06_spos_p1 = (0.039 + (3.5 + fX));  //Christoph
	ss06_sbw_p1 = 2;
	ss06_sarea_p1 = 0;
	ss06_seta_p1 = 0;
	  
	ss06_ke_p1 = 50;
//	ss06_kpos_p1 = (0.00798 + (2 + fY));//detectors couldn't be shifted in simulation, they are shifted here 0.00798
        ss06_kpos_p1 = (2 + fY);
//	ss06_kpos_p1 = (0.01035 + (2 + fY));  //Christoph
	ss06_kbw_p1 = 2;
	ss06_karea_p1 = 0;
	ss06_keta_p1 = 0;
//	cout<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID<<endl;   
      }
      
      if (DetID==16 && l>1)
      {
        ss06_se_p2 = 50;
//	ss06_spos_p2 = (0.038495 + (3.5 + fX));//detectors couldn't be shifted in simulation, they are shifted here -0.038495
	ss06_spos_p2 = (3.5 + fX);
//	ss06_spos_p1 = (0.039 + (3.5 + fX));  //Christoph
	ss06_sbw_p2 = 2;
	ss06_sarea_p2 = 0;
	ss06_seta_p2 = 0;
	  
	ss06_ke_p2 = 50;
//	ss06_kpos_p2 = (0.00798 + (2 + fY));//detectors couldn't be shifted in simulation, they are shifted here 0.00798
	ss06_kpos_p2 = (2 + fY);
//	ss06_kpos_p2 = (0.01035 + (2 + fY));  //Christoph
	ss06_kbw_p2 = 2;
	ss06_karea_p2 = 0;
	ss06_keta_p2 = 0;
//	cout<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID<<endl;   
      }
      
    }     

 }
   
  


AddHit(ss03_se_p1,ss03_spos_p1,ss03_sbw_p1,ss03_sarea_p1,
ss03_seta_p1,ss03_ke_p1,ss03_kpos_p1,ss03_kbw_p1,ss03_karea_p1,ss03_keta_p1,ss06_se_p1,ss06_spos_p1,ss06_sbw_p1,
ss06_sarea_p1,ss06_seta_p1,ss06_ke_p1,ss06_kpos_p1,ss06_kbw_p1,ss06_karea_p1,ss06_keta_p1,ss03_se_p2,ss03_spos_p2,
ss03_sbw_p2,ss03_sarea_p2,ss03_seta_p2,ss03_ke_p2,ss03_kpos_p2,ss03_kbw_p2,ss03_karea_p2,ss03_keta_p2,ss06_se_p2,
ss06_spos_p2,ss06_sbw_p2,ss06_sarea_p2,ss06_seta_p2,ss06_ke_p2,ss06_kpos_p2,ss06_kbw_p2,ss06_karea_p2,ss06_keta_p2);


}
// -------------------------------------------------------------------------

void R3BTra2pDigitizer::Reset(){
// Clear the structure
//   cout << " -I- Digit Reset() called " << endl;

   
 if (fTra2pDigi ) fTra2pDigi->Clear();

}   

void R3BTra2pDigitizer::Finish()
{
// Write control histograms

   
}

R3BTra2pDigi* R3BTra2pDigitizer::AddHit(Double_t ss03_se_p1,
Double_t ss03_spos_p1,Double_t ss03_sbw_p1,Double_t ss03_sarea_p1,Double_t ss03_seta_p1,Double_t ss03_ke_p1,
Double_t ss03_kpos_p1,Double_t ss03_kbw_p1,Double_t ss03_karea_p1,Double_t ss03_keta_p1,Double_t ss06_se_p1,
Double_t ss06_spos_p1,Double_t ss06_sbw_p1,Double_t ss06_sarea_p1,Double_t ss06_seta_p1,Double_t ss06_ke_p1,
Double_t ss06_kpos_p1,Double_t ss06_kbw_p1,Double_t ss06_karea_p1,Double_t ss06_keta_p1,Double_t ss03_se_p2,
Double_t ss03_spos_p2,Double_t ss03_sbw_p2,Double_t ss03_sarea_p2,Double_t ss03_seta_p2,Double_t ss03_ke_p2,
Double_t ss03_kpos_p2,Double_t ss03_kbw_p2,Double_t ss03_karea_p2,Double_t ss03_keta_p2,Double_t ss06_se_p2,
Double_t ss06_spos_p2,Double_t ss06_sbw_p2,Double_t ss06_sarea_p2,Double_t ss06_seta_p2,Double_t ss06_ke_p2,
Double_t ss06_kpos_p2,Double_t ss06_kbw_p2,Double_t ss06_karea_p2,Double_t ss06_keta_p2){   
  TClonesArray& clref = *fTra2pDigi;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) R3BTra2pDigi(ss03_se_p1,ss03_spos_p1,ss03_sbw_p1,ss03_sarea_p1,
ss03_seta_p1,ss03_ke_p1,ss03_kpos_p1,ss03_kbw_p1,ss03_karea_p1,ss03_keta_p1,ss06_se_p1,ss06_spos_p1,ss06_sbw_p1,
ss06_sarea_p1,ss06_seta_p1,ss06_ke_p1,ss06_kpos_p1,ss06_kbw_p1,ss06_karea_p1,ss06_keta_p1,ss03_se_p2,ss03_spos_p2,
ss03_sbw_p2,ss03_sarea_p2,ss03_seta_p2,ss03_ke_p2,ss03_kpos_p2,ss03_kbw_p2,ss03_karea_p2,ss03_keta_p2,ss06_se_p2,
ss06_spos_p2,ss06_sbw_p2,ss06_sarea_p2,ss06_seta_p2,ss06_ke_p2,ss06_kpos_p2,ss06_kbw_p2,ss06_karea_p2,ss06_keta_p2);
 
}



//R3BTraDigi* R3BTra2pDigitizer::AddHit(
 // return new(clref[size]) R3BTra2pDigi();
//}

ClassImp(R3BTra2pDigitizer)