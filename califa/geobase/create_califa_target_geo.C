/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// FINAL CALIFA CARREL + iPHOS
// Usage after loading:
// create_califa_geo_selector(A, B)
// A: known alignment for some experiments
//    [nominal, s522, s509, s455, s515, s494, s444, s467]
// B: geoTag. Name tagging the output root file.
// Requires file CLF-ALL-oneCrystal.txt for crystals coordinates
// Requires file CLF-ALL-onePart.txt for alveoli coordinates

#include "TGeoManager.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>

const int nbcrystals = 2544;

Bool_t isCrystalInstalled(Int_t alvType, Int_t alveolusCopy, Int_t instCry[]);

// Create Matrix Unity
TGeoRotation *fGlobalRot = new TGeoRotation();

// Create a null translation
TGeoTranslation *fGlobalTrans = new TGeoTranslation();

//TGeoVolume* gTop;
TGeoRotation* fRefRot = NULL;
TGeoManager* gGeoMan = NULL;
TGeoVolume* pWorld = NULL;


Double_t fThetaX = 0.;
Double_t fThetaY = 0.;
Double_t fThetaZ = 0.;
Double_t fPhi   = 0.;
Double_t fTheta = 0.;
Double_t fPsi   = 0.;
Double_t fX = 0.;
Double_t fY = 0.;
Double_t fZ = 0.;
Bool_t fLocalTrans = kFALSE;
Bool_t fLabTrans = kFALSE;

TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef);
void ConstructGeometry(TGeoMedium *pMedVac, TGeoMedium *pMed1,
		       TGeoMedium *pMed2, TGeoMedium *pMed3,
		       TGeoMedium *pMed4, TGeoMedium *pMed5,
		       TGeoMedium *pMedSi, TGeoMedium *pMedGold, TGeoMedium *pMedPolyethylene, TGeoMedium *pMed41, TGeoMedium *pMed42, TGeoMedium *pMed43);



void create_califa_geo_selector(const char* expNumber = "nominal", TString geoTag = "full")
{

    // -------   Load media from media file   -------------------------
    FairGeoLoader* geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    TString geoPath = gSystem->Getenv("VMCWORKDIR");
    TString medFile = geoPath + "/geometry/media_r3b.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    gGeoMan = gGeoManager;

    TString installedCrystalsFile = "califa_AllCrystalsInstalled.txt";
    if ((strncmp(expNumber, "s091", 4) == 0) || (strncmp(expNumber, "s118", 4) == 0))
    {
        installedCrystalsFile = "califa_InstalledCrystals_Jan2024.txt";
        geoTag = "v2024.1";
    }
    else if ((strncmp(expNumber, "s522", 4) == 0) || (strncmp(expNumber, "s515", 4) == 0) ||
             (strncmp(expNumber, "s509", 4) == 0) || (strncmp(expNumber, "s494", 4) == 0) ||
             (strncmp(expNumber, "s455", 4) == 0) || (strncmp(expNumber, "s444", 4) == 0))
    {
        installedCrystalsFile = "califa_InstalledCrystals_March2021.txt";
        geoTag = "v2021.3";
    }
    else if (strncmp(expNumber, "s467", 4) == 0)
    {
        installedCrystalsFile = "califa_InstalledCrystals_Nov2019.txt";
        geoTag = "v2019.11";
    }
    else if (strncmp(expNumber, "nominal", 7) == 0)
    {
        // Nominal configuration, nothing to say.
        std::cout << std::endl;
        std::cout << "\033[35m"
                  << "Loading the ideal nominal configuration."
                  << " \033[0m" << std::endl;
    }
    else
    {
        std::cout << std::endl
                  << "\033[35m WARNING: Using a non-standard experiment number (" << expNumber
                  << ") therefore the CALIFA geometry will be generated with all crystals using the file "
                     "<califa_AllCrystalsInstalled.txt> \033[0m"
                  << std::endl;
        std::cout << "\033[33m"
                  << "  -> Is this really what you want??"
                  << " \033[0m" << std::endl;
    }

    /*
        if ((strncmp(expNumber, "s522", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s509", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s455", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s515", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s494", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s444", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_March2021.txt", 38) != 0) ||
            (strncmp(expNumber, "s467", 4) == 0 &&
             strncmp(installedCrystalsFile, "califa_InstalledCrystals_Nov2019.txt", 38) != 0))
        {
            std::cout << std::endl
                      << "WARNING: Using a non-standard combination of experiment number (" << expNumber
                      << ") and installed crystals file (" << installedCrystalsFile << ")" << std::endl;
            std::cout << "\033[33m"
                      << "  -> Is this really what you want??"
                      << " \033[0m" << std::endl;
        }*/

    TString fFilePath = geoPath + "/califa/geobase/files/" + installedCrystalsFile;
    fFilePath.ReplaceAll("./", "/");
    fFilePath.ReplaceAll("//", "/");

    std::cout << std::endl
              << " geoTag: "
              << "\033[33m" << geoTag << " \033[0m" << std::endl;
    std::cout << " alignment (nominal or experiment code): "
              << "\033[33m" << expNumber << " \033[0m" << std::endl;
    std::cout << " installed Crystals File: "
              << "\033[33m" << fFilePath << " \033[0m" << std::endl
              << std::endl;

    if (strncmp(expNumber, "s091", 4) == 0 || strncmp(expNumber, "s118", 4) == 0 ||
        strncmp(expNumber, "s522", 4) == 0 || strncmp(expNumber, "s509", 4) == 0 ||
        strncmp(expNumber, "s494", 4) == 0 || strncmp(expNumber, "s444", 4) == 0 || strncmp(expNumber, "s467", 4) == 0)
    {
        std::cout << "\033[31m No alignment data for this experiment yet \033[0m " << std::endl << std::endl;
    }

    ifstream wc1, in1, in2, in3;
    wc1.open(fFilePath.Data());
    in1.open((geoPath + "/califa/geobase/files/CLF-ALL-onePart.txt").Data());
    in2.open((geoPath + "/califa/geobase/files/CLF-ALL-oneCrystal.txt").Data());
    in3.open((geoPath + "/califa/geobase/files/CLF-CEPA-CSI_Crystal.txt").Data());

    if (wc1.fail())
    {
        std::cout << "\033[31m FATAL ERROR: Unable to open file or read problem in \033[0m "
                  << "\033[31m" << installedCrystalsFile << " \033[0m" << std::endl
                  << std::endl;
        return;
    }

    Int_t installedCrystals[nbcrystals];
    Int_t read = 0;
    Int_t crycounter = 0;
    for (Int_t i = 0; i < nbcrystals; i++)
        installedCrystals[i] = 0;
    while (1)
    { // reading the file with all alveoli vertices
        wc1 >> read;
        if (!wc1.good())
        {
            break;
        }
        if (crycounter > (nbcrystals - 1))
        {
            std::cout << "\033[31m FATAL ERROR: more lines than crystals in the file \033[0m "
                      << "\033[31m" << installedCrystalsFile << " \033[0m" << std::endl
                      << std::endl;
            return;
        }
        installedCrystals[crycounter] = read;
        crycounter++;
    }

    // -------   Geometry file name (output)   ------------------------
    TString geoFileName = geoPath + "/geometry/califa_";
    geoFileName = geoFileName + geoTag + ".geo.root";

    // -----------------   Get and create the required media    -------
    FairGeoMedia* geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    FairGeoMedium* mAir = geoMedia->getMedium("Air");
    if (!mAir)
        Fatal("Main", "FairMedium Air not found");
    geoBuild->createMedium(mAir);
    TGeoMedium* pAirMedium = gGeoMan->GetMedium("Air");
    if (!pAirMedium)
        Fatal("Main", "Medium Air not found");

    FairGeoMedium* mCsI = geoMedia->getMedium("CsI");
    if (!mCsI)
        Fatal("Main", "FairMedium CsI not found");
    geoBuild->createMedium(mCsI);
    TGeoMedium* pCsIMedium = gGeoMan->GetMedium("CsI");
    if (!pCsIMedium)
        Fatal("Main", "Medium CsI not found");

    FairGeoMedium* mCar = geoMedia->getMedium("CarbonFibre");
    if (!mCar)
        Fatal("Main", "FairMedium CarbonFibre not found");
    geoBuild->createMedium(mCar);
    TGeoMedium* pCarbonFibreMedium = gGeoMan->GetMedium("CarbonFibre");
    if (!pCarbonFibreMedium)
        Fatal("Main", "Medium CarbonFibre not found");

    FairGeoMedium* mTfl = geoMedia->getMedium("Tefflon");
    if (!mTfl)
        Fatal("Main", "FairMedium Tefflon not found");
    geoBuild->createMedium(mTfl);
    TGeoMedium* pWrappingMedium = gGeoMan->GetMedium("Tefflon");
    if (!pWrappingMedium)
        Fatal("Main", "Medium Tefflon not found");

  /*------------------------------------*/ //FOOT
    FairGeoMedium* mVac = geoMedia->getMedium("vacuum");
  if (!mVac)
    Fatal("Main", "FairMedium vacuum not found");
  geoBuild->createMedium(mVac);
  TGeoMedium* pMedVac = gGeoMan->GetMedium("vacuum");
  if (!pMedVac)
    Fatal("Main", "Medium vacuum not found");

  FairGeoMedium* mCu = geoMedia->getMedium("copper");
  if (!mCu)
    Fatal("Main", "FairMedium copper not found");
  geoBuild->createMedium(mCu);
  TGeoMedium* pMed1 = gGeoMan->GetMedium("copper");
  if (!pMed1)
    Fatal("Main", "Medium copper not found");

  FairGeoMedium* mAl = geoMedia->getMedium("aluminium");
  if (!mAl)
    Fatal("Main", "FairMedium aluminium not found");
  geoBuild->createMedium(mAl);
  TGeoMedium* pMed2 = gGeoMan->GetMedium("aluminium");
  if (!pMed2)
    Fatal("Main", "Medium aluminium not found");

  FairGeoMedium* mSteel = geoMedia->getMedium("Steel");
  if (!mSteel)
    Fatal("Main", "FairMedium Steel not found");
  geoBuild->createMedium(mSteel);
  TGeoMedium* pMed3 = gGeoMan->GetMedium("Steel");
  if (!pMed3)
    Fatal("Main", "Medium Steel not found");

  FairGeoMedium* mMylar = geoMedia->getMedium("mylar");
  if (!mMylar)
    Fatal("Main", "FairMedium mylar not found");
  geoBuild->createMedium(mMylar);
  TGeoMedium* pMed4 = gGeoMan->GetMedium("mylar");
  if (!pMed4)
    Fatal("Main", "Medium mylar not found");

  FairGeoMedium* mH2 = geoMedia->getMedium("H2");
  if (!mH2)
    Fatal("Main", "FairMedium H2 not found");
  geoBuild->createMedium(mH2);
  TGeoMedium* pMed5 = gGeoMan->GetMedium("H2");
  if (!pMed5)
    Fatal("Main", "Medium H2 not found");

  FairGeoMedium* mSi = geoMedia->getMedium("silicon");
  if (!mSi)
    Fatal("Main", "FairMedium silicon not found");
  geoBuild->createMedium(mSi);
  TGeoMedium* pMedSi = gGeoMan->GetMedium("silicon");
  if (!pMedSi)
    Fatal("Main", "Medium silicon not found");

  FairGeoMedium* mGold      = geoMedia->getMedium("gold");
  if ( ! mGold ) Fatal("Main", "FairMedium gold not found");
  geoBuild->createMedium(mGold);
  TGeoMedium* pMedGold = gGeoMan->GetMedium("gold");
  if ( ! pMedGold ) Fatal("Main", "Medium gold not found");

  FairGeoMedium* mPolyethylene      = geoMedia->getMedium("polyethylene");
  if ( ! mPolyethylene ) Fatal("Main", "FairMedium polyethylene not found");
  geoBuild->createMedium(mPolyethylene);
  TGeoMedium* pMedPolyethylene = gGeoMan->GetMedium("polyethylene");
  if ( ! pMedPolyethylene ) Fatal("Main", "Medium polyethylene not found");

  FairGeoMedium *mCCooling = geoMedia->getMedium("carboncooling");
  if (!mCCooling)
    Fatal("Main", "FairMedium carboncooling not found");
  geoBuild->createMedium(mCCooling);
  TGeoMedium *pMed42 = gGeoMan->GetMedium("carboncooling");
  if (!pMed42)
    Fatal("Main", "Medium carboncooling not found");
  
  FairGeoMedium *mPolyamide = geoMedia->getMedium("polyamide");
  if (!mPolyamide)
    Fatal("Main", "FairMedium polyamide not found");
  geoBuild->createMedium(mPolyamide);
  TGeoMedium *pMed43 = gGeoMan->GetMedium("polyamide");
  if (!pMed43)
    Fatal("Main", "Medium polyamide not found");
    
  FairGeoMedium *mKapton = geoMedia->getMedium("kapton");
  if (!mKapton)
    Fatal("Main", "FairMedium kapton not found");
  geoBuild->createMedium(mKapton);
  TGeoMedium *pMed41 = gGeoMan->GetMedium("kapton");
  if (!pMed41)
    Fatal("Main", "Medium kapton not found");  
  
  
    // ----------------------------------------------------------------

    // --------------   Create geometry and top volume  ---------------
    gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    gGeoMan->SetName("CALIFAgeom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    gGeoMan->SetTopVolume(top);
    // ----------------------------------------------------------------

    // Defintion of the Mother Volume
    auto pCBWorldOut = new TGeoTube("Califa_boxOut",
                                    0.,               // Rmin
                                    59.,              // Rmax
                                    //(55. + 70) / 2.); // half length
                                    (55. + 170)); // half length
    auto trans0 = new TGeoCombiTrans("trans0", 0., 0., 7.5-100, fRefRot);
    trans0->RegisterYourself();

    TGeoShape* pCBWorldIn1 = new TGeoTube("Califa_Centerpart1", // hole to accommodate the tracker
                                          0.,                   // Rmin
                                          27.,                  // Rmax
                                          //130 / 2.);            // half length
                                          130-100);            // half length
    auto trans1 = new TGeoCombiTrans("trans1", 0., 0., 7.5+100, fRefRot);
    trans1->RegisterYourself();

    auto Part2 = new TGeoSphere("Part2", 38., 73.5, 7., 90., 0., 360.);
    auto trans2 = new TGeoCombiTrans("trans2", 0, 0, 0., fRefRot);
    trans2->RegisterYourself();

    auto pCBWorld = new TGeoCompositeShape("Califa_box", "Califa_boxOut:trans0-Califa_Centerpart1:trans1+Part2:trans2");

    pWorld = new TGeoVolume("CalifaWorld", pCBWorld, pAirMedium);
    auto t0 = new TGeoCombiTrans(0., 0., 0., fRefRot);
    top->AddNode(pWorld, 0, t0);

    // FINAL CALIFA CARREL + iPHOS VERSION (SINCE NOV 2019)
    const Int_t N_ALV_TYPES = 23; // alveolar structures
    const Int_t N_CRY_TYPES = 85; // crystal elements
    // CALIFA CEPA USC VERSION (JAN 2024) PARAMETERS
    const Int_t N_ALV_TYPES_CEPA = 3;  // alveolar structures
    const Int_t N_CRY_TYPES_CEPA = 14; // crystal elements

    Int_t counter = 0;
    Float_t x, y, z;

    Double_t wrapping_thickness = 0.0065;            // in cm.
    const Double_t wrapping_thickness_CEPA = 0.0080; // (80 microns recomended by E. Casarejos)
    const Double_t cf_thickness_CEPA = 0.0300;       // carbon fiber (300 microns recomended by E. Casarejos)
    const Double_t security_margin =
        0.0110; // between wrapping and carbon fiber inner volume (minimum 40 microns, depends on crystal_reduction)
    const Double_t crystal_reduction = -0.0330; // minimum -220 microns as there was a small overlap among sectors
    // the last two parameters depends one on the other... the minimum values recommended in the comment
    // avoids a given type of extrusion (security_margin) or overlap (crystal_reduction), but to avoid simultaneously
    // both problems below 1 micron, the minimum settings are security_margin = 0.0110 and crystal_reduction = -0.0330

    // target reference in mm. OFFSET INFO FROM UVIGO
    TVector3 target_ref(4.1, 2304.0809, 325.0);
    TVector3 target_ref_CEPA(0.0, 0.0, 24.9);

    // 23 geometries, 8 vertices, outer and inner: (23*8*2)
    TVector3 points[N_ALV_TYPES * 8 * 2];
    TVector3 points_local[N_ALV_TYPES * 8 * 2];

    // 18 alv with 4 cry, 4 with 3 cry, 1 with 1 cry, 8 vertices each: (15*4+4*3+3*4+1)*8=85*8
    TVector3 points_cry[N_CRY_TYPES * 8];
    TVector3 points_cry_local[N_CRY_TYPES * 8];

    // CEPA 3 geometries, 8 vertices = 24
    TVector3 points_CEPA[N_ALV_TYPES_CEPA * 8];
    TVector3 points_local_CEPA[N_ALV_TYPES_CEPA * 8];
    TVector3 points_inn_CEPA[N_ALV_TYPES_CEPA * 8];
    TVector3 points_inn_local_CEPA[N_ALV_TYPES_CEPA * 8];

    // 14 geometries, 8 vertices = 112
    TVector3 points_cry_CEPA[N_CRY_TYPES_CEPA * 8];
    TVector3 points_cry_wrap_CEPA[N_CRY_TYPES_CEPA * 8];
    TVector3 points_cry_local_CEPA[N_CRY_TYPES_CEPA * 8];
    TVector3 points_cry_wrap_local_CEPA[N_CRY_TYPES_CEPA * 8];

    while (1)
    { // reading the file with all alveoli vertices
        in1 >> x >> y >> z;
        if (!in1.good())
            break;
        points[counter].SetXYZ(
            (x - target_ref.X()) / 10, (y - target_ref.Y()) / 10, (z - target_ref.Z()) / 10); // in cm;
        if (points[counter].X() > 100 || points[counter].Y() > 100 || points[counter].Z() > 100)
            cout << "warn: points exceed top volume!!" << endl;
        counter++;
        // printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    }
    if (counter != N_ALV_TYPES * 8 * 2)
        cout << "PROBLEM! Counter1=" << counter << endl;
    counter = 0;

    while (1)
    { // reading the file with all crystal vertices
        in2 >> x >> y >> z;
        if (!in2.good())
            break;
        points_cry[counter].SetXYZ(
            (x - target_ref.X()) / 10, (y - target_ref.Y()) / 10, (z - target_ref.Z()) / 10); // in cm;
        if (points_cry[counter].X() > 100 || points_cry[counter].Y() > 100 || points_cry[counter].Z() > 100)
            cout << "warn: points exceed top volume!!" << endl;
        counter++;
        // printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    }
    if (counter != N_CRY_TYPES * 8)
        cout << "PROBLEM! Counter2=" << counter << endl;
    counter = 0;

    while (1)
    { // reading the file with CEPA vertices
        in3 >> x >> y >> z;
        if (!in3.good())
            break;
        points_cry_CEPA[counter].SetXYZ(
            (x - target_ref_CEPA.X()) / 10, (y - target_ref_CEPA.Y()) / 10, (z - target_ref_CEPA.Z()) / 10); // in cm;
        if (points_cry_CEPA[counter].X() > 100 || points_cry_CEPA[counter].Y() > 100 ||
            points_cry_CEPA[counter].Z() > 100)
            cout << "WARNING: points exceed top volume!!" << endl;
        counter++;
        // printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    }
    if (counter != 112)
        cout << "PROBLEM! Counter3=" << counter << endl;
    counter = 0;

    // BARREL+iPHOS PART
    //  centers of faces
    TVector3 center[N_ALV_TYPES * 2 * 2]; // 23 geometries, 2 face centers, outer and inner (23*2*2)
    TVector3 x_uni[N_ALV_TYPES * 2 * 2];  // unit vectors for each face
    TVector3 y_uni[N_ALV_TYPES * 2 * 2];
    TVector3 z_uni[N_ALV_TYPES * 2 * 2];
    TVector3 center_cry[N_CRY_TYPES * 2]; // 85 types of crystals, 2 face centers (85*2)
    TVector3 x_uni_cry[N_CRY_TYPES * 2];  // unit vectors for each face for crystals
    TVector3 y_uni_cry[N_CRY_TYPES * 2];
    TVector3 z_uni_cry[N_CRY_TYPES * 2];
    TRotation rot[N_ALV_TYPES * 2 * 2 + 1]; // calculated in each face, but only 23 are really different if all is ok
    // NOTE: the +1 at the end solves a problem in some linux systems where there is a overwritting of some values of
    // this matrix. We could not identify the origin of this (array out of bounds??) which is only present in some
    // systems.
    TRotation rot_cry[N_CRY_TYPES * 2]; // only a few are really different if all is ok
    // volume centers
    TVector3 alv_cm[N_ALV_TYPES * 2];     // 23 geometries, outer and inner
    TVector3 alv_cm_rot[N_ALV_TYPES * 2]; // 23 geometries, outer and inner after final rotation
    TVector3 cry_cm[N_CRY_TYPES];         // 85 types of crystals

    // The center of the faces are first calculated. Then, the unit vectors defining the axis in each faces
    // Third, the rotation moving from the lab system to the unit vectors previously found. To define the
    // volume in Arb8 style, we need the 8 corners in the local frustrum coordinates. Then, we should express
    // the vertices in the coordinate system of the volume center of mass (cm)
    for (Int_t i = 0; i < N_ALV_TYPES * 2 * 2; i++)
    { // for 23 geometries, 2 face centers, outer and inner (23*2*2)
        center[i] = points[i * 4] + points[i * 4 + 1] + points[i * 4 + 2] + points[i * 4 + 3]; // face centers
        center[i] *= 0.25;                                                                     // face centers
        // center[i].Print();
        z_uni[i] = (points[i * 4 + 1] - points[i * 4]).Cross(points[i * 4 + 2] - points[i * 4 + 1]);
        z_uni[i] = z_uni[i].Unit(); // normal to face center
        x_uni[i] = points[i * 4 + 1] - points[i * 4];
        if (i > 75)
        { // for the irregular alveoli at the end of the iPhos, the X axis is taken from other vertices
            x_uni[i] = points[i * 4 + 2] - points[i * 4 + 3];
        }
        x_uni[i] = x_uni[i].Unit();          // unit along X
        y_uni[i] = z_uni[i].Cross(x_uni[i]); // unit along Y

        // x_uni[i].Print();  y_uni[i].Print();   z_uni[i].Print();

        // calculate rotation matrix for the 23 geometries (should be repeated 4 times, just checking)
        rot[i].SetZAxis(z_uni[i], x_uni[i]);
        // cout << "rot["<< i<< "]" << rot[i].XX() << " "<< rot[i].XY() << " "<< rot[i].XZ()<< endl;
        // cout << "rot["<< i<< "]" << rot[i].YX() << " "<< rot[i].YY() << " "<< rot[i].YZ()<< endl;
        // cout << "rot["<< i<< "]" << rot[i].ZX() << " "<< rot[i].ZY() << " "<< rot[i].ZZ()<< endl;

        // calculate vertices coordinates in local alveoli coordinates (not valid for the definition of the Arb8, see
        // below) points_local[i*4]=rot[i].Inverse()*(points[i*4]-center[i]);
        // points_local[i*4+1]=rot[i].Inverse()*(points[i*4+1]-center[i]);
        // points_local[i*4+2]=rot[i].Inverse()*(points[i*4+2]-center[i]);
        // points_local[i*4+3]=rot[i].Inverse()*(points[i*4+3]-center[i]);
        // cout<< "Points in face coordinates: "<< endl;
        // points_local[i*4].Print();points_local[i*4+1].Print();points_local[i*4+2].Print();points_local[i*4+3].Print();
    }
    for (Int_t i = 0; i < N_ALV_TYPES * 2; i++)
    { // for 23 geometries, outer and inner (23*2)
        alv_cm[i] = center[2 * i] + center[2 * i + 1];
        alv_cm[i] *= 0.5; // volume center for all alv (outer and inner)
        // alv_cm[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each alveolus
            points_local[i * 8 + j] = rot[2 * i].Inverse() * (points[i * 8 + j] - alv_cm[i]);
            // cout<< "Points in cm coordinates: "<< endl;
            // points_local[i*8+j].Print();
        }
    }

    // same procedure for crystals
    for (Int_t i = 0; i < N_CRY_TYPES * 2; i++)
    { // 85 types of crystals, 2 face centers (85*2)
        center_cry[i] =
            points_cry[i * 4] + points_cry[i * 4 + 1] + points_cry[i * 4 + 2] + points_cry[i * 4 + 3]; // face centers
        center_cry[i] *= 0.25;                                                                         // face centers
        // center_cry[i].Print();
        z_uni_cry[i] = (points_cry[i * 4 + 1] - points_cry[i * 4]).Cross(points_cry[i * 4 + 2] - points_cry[i * 4 + 1]);
        z_uni_cry[i] = z_uni_cry[i].Unit(); // normal to face center
        x_uni_cry[i] = points_cry[i * 4 + 1] - points_cry[i * 4];
        if (i > 145)
        { // for the irregular crystals at the end of the iPhos, the X axis is taken from other vertices
            x_uni[i] = points_cry[i * 4 + 2] - points_cry[i * 4 + 3];
        }
        x_uni_cry[i] = x_uni_cry[i].Unit();              // unit along X
        y_uni_cry[i] = z_uni_cry[i].Cross(x_uni_cry[i]); // unit along Y
        // x_uni_cry[i].Print();  y_uni_cry[i].Print();    z_uni_cry[i].Print();

        // calculate rotation matrix for the 85 types of crystals (should be repeated 4 times, just checking)
        rot_cry[i].SetZAxis(z_uni_cry[i], x_uni_cry[i]);

        // calculate vertices coordinates in local crystal coordinates (not valid for the definition of the Arb8, see
        // below) points_cry_local[i*4]=rot_cry[i].Inverse()*(points_cry[i*4]-center_cry[i]);
        // points_cry_local[i*4+1]=rot_cry[i].Inverse()*(points_cry[i*4+1]-center_cry[i]);
        // points_cry_local[i*4+2]=rot_cry[i].Inverse()*(points_cry[i*4+2]-center_cry[i]);
        // points_cry_local[i*4+3]=rot_cry[i].Inverse()*(points_cry[i*4+3]-center_cry[i]);
        // cout<< "Points in face coordinates: "<< endl;
        // points_cry_local[i*4].Print();points_cry_local[i*4+1].Print();
        // points_cry_local[i*4+2].Print();points_cry_local[i*4+3].Print();
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    { // 85 types of crystals
        cry_cm[i] = center_cry[2 * i] + center_cry[2 * i + 1];
        cry_cm[i] *= 0.5; // volume center for each cry
        // cry_cm[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each crystal
            points_cry_local[i * 8 + j] = rot_cry[2 * i].Inverse() * (points_cry[i * 8 + j] - cry_cm[i]);
            // reducing 1mm the first crystal vertices to avoid collision with inner volume,
            // while the first crystal definition is not fixed
            if (i == 0)
            {
                if (points_cry_local[i * 8 + j].X() > 0)
                {
                    points_cry_local[i * 8 + j].SetX(points_cry_local[i * 8 + j].X() - 0.1);
                }
                else
                {
                    points_cry_local[i * 8 + j].SetX(points_cry_local[i * 8 + j].X() + 0.1);
                }
                if (points_cry_local[i * 8 + j].Y() > 0)
                {
                    points_cry_local[i * 8 + j].SetY(points_cry_local[i * 8 + j].Y() - 0.1);
                }
                else
                {
                    points_cry_local[i * 8 + j].SetY(points_cry_local[i * 8 + j].Y() + 0.1);
                }
            }
            if (i > 72)
            { // for the irregular crystals at the end of the iPhos, rotation is taken from the corresponding alveolus
                points_cry_local[i * 8 + j] =
                    rot[4 * (Int_t)((i - 73) / 3) + 78].Inverse() * (points_cry[i * 8 + j] - cry_cm[i]);
            }
            // cout<< "Points in cm coordinates.  crystal " << i << ",  vertex "<< j<< endl;
            // points_cry_local[i*8+j].Print();
        }
    }

    // location of the crystals in the alveoli
    TVector3 cry_position[N_CRY_TYPES];       // 85 types of crystals
    TVector3 cry_position_local[N_CRY_TYPES]; // 85 types of crystals

    cry_position[0] = cry_cm[0] - alv_cm[1]; // first alveolus with a single crystal wrt inner alv
    cry_position_local[0] = rot[1].Inverse() * (cry_cm[0] - alv_cm[0]);
    // Relative Crystal rotation in each alveoli. Obtained from the crystal unit vector in
    for (Int_t i = 1; i < N_CRY_TYPES; i++)
    {
        if (i < 73)
        { // four crystals per alv
            // cout << i <<  " " << 2*(Int_t)((i-1)/4)+3 <<  " "<< 4*(Int_t)((i-1)/4)+6 << endl;
            cry_position[i] = cry_cm[i] - alv_cm[2 * (Int_t)((i - 1) / 4) + 3];
            // cry_position[i].Print();

            // cout << "cry_cm ";cry_cm[i].Print();
            // cout << "alv_cm ";alv_cm[(Int_t)((i-1)/4)+1].Print();
            // cout << "rot["<<4 * (Int_t)((i - 1) / 4) + 6<<"]" << rot[4 * (Int_t)((i - 1) / 4) + 6].XX() << " "<<
            // rot[4 * (Int_t)((i - 1) / 4) + 6].XY() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].XZ()<< endl; cout <<
            // "rot["<<4 * (Int_t)((i - 1) / 4) + 6<<"]" << rot[4 * (Int_t)((i - 1) / 4) + 6].YX() << " "<< rot[4 *
            // (Int_t)((i - 1) / 4) + 6].YY() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].YZ()<< endl; cout << "rot["<<4
            // * (Int_t)((i - 1) / 4) + 6<<"]" << rot[4 * (Int_t)((i - 1) / 4) + 6].ZX() << " "<< rot[4 * (Int_t)((i -
            // 1) / 4) + 6].ZY() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].ZZ()<< endl; cout << "rot_inverse " << rot[4
            // * (Int_t)((i - 1) / 4) + 6].Inverse().XX() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].Inverse().XY()<< "
            // "<<rot[4 * (Int_t)((i - 1) / 4) + 6].Inverse().XZ()<< endl; cout << "rot_inverse " << rot[4 * (Int_t)((i
            // - 1) / 4) + 6].Inverse().YX() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].Inverse().YY()<< " "<<rot[4 *
            // (Int_t)((i - 1) / 4) + 6].Inverse().YZ()<< endl; cout << "rot_inverse " << rot[4 * (Int_t)((i - 1) / 4) +
            // 6].Inverse().ZX() << " "<< rot[4 * (Int_t)((i - 1) / 4) + 6].Inverse().ZY()<< " "<<rot[4 * (Int_t)((i -
            // 1) / 4) + 6].Inverse().ZZ()<< endl;

            cry_position_local[i] =
                rot[4 * (Int_t)((i - 1) / 4) + 6].Inverse() * (cry_cm[i] - alv_cm[2 * (Int_t)((i - 1) / 4) + 3]);
            // TODO, KNOWN ISSUE TO SOLVE! Moving up 15mm all crystals in BARREL and iPhos to avoid collisions!!!
            // All collisions were below 200 microns, but it should be checked in the original CAD files
            cry_position_local[i].SetZ(cry_position_local[i].Z() + 1.5);
            // cout << "local " ;cry_position_local[i].Print();
        }
        else
        { // last four alveoli with three crystals per alv
            // cout << i <<  " " << 2*(Int_t)((i-73)/3)+39 <<  " " << 4*(Int_t)((i-73)/3)+78 << endl;
            cry_position[i] = alv_cm[2 * (Int_t)((i - 73) / 3) + 39] - cry_cm[i];
            cry_position_local[i] =
                rot[4 * (Int_t)((i - 73) / 3) + 78].Inverse() * (cry_cm[i] - alv_cm[2 * (Int_t)((i - 73) / 3) + 39]);
            // TODO, KNOWN ISSUE TO SOLVE! Moving up 15mm all crystals in BARREL and iPhos to avoid collisions!!!
            // There is no collision for the alveoli with 3 crystals... anyway, moving them for coherence
            cry_position_local[i].SetZ(cry_position_local[i].Z() + 1.5);
            // cry_position[i].Print();
            // cry_position_local[i].Print();
        }
    }

    // Redefinition of vertices for the construction of the Alveoli, using TGeoArb8
    Double_t* vertices_Alv[N_ALV_TYPES]; // 23 geometries
    Double_t* vertices_inner_Alv[N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        vertices_Alv[i] = new Double_t[16];
        vertices_inner_Alv[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
            if ((4 - j) > 0)
            { // reversing order for being clockwise filling TGeoArb8
                vertices_Alv[i][2 * j] = points_local[16 * i + 3 - j].X();
                vertices_Alv[i][2 * j + 1] = points_local[16 * i + 3 - j].Y();
                vertices_inner_Alv[i][2 * j] = points_local[16 * i + 3 - j + 8].X();
                vertices_inner_Alv[i][2 * j + 1] = points_local[16 * i + 3 - j + 8].Y();
            }
            else
            {
                vertices_Alv[i][2 * j] = points_local[16 * i + 11 - j].X();
                vertices_Alv[i][2 * j + 1] = points_local[16 * i + 11 - j].Y();
                vertices_inner_Alv[i][2 * j] = points_local[16 * i + 11 - j + 8].X();
                vertices_inner_Alv[i][2 * j + 1] = points_local[16 * i + 11 - j + 8].Y();
            }
        }
    }

    // Redefinition of vertices for the construction of the Crystals, using TGeoArb8
    Double_t* vertices_Cry[N_CRY_TYPES];      // 680/8=85 (85 crystal types)
    Double_t* vertices_Cry_Wrap[N_CRY_TYPES]; // 680/8=85 (85 crystal types)
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        vertices_Cry[i] = new Double_t[16];
        vertices_Cry_Wrap[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
            if ((4 - j) > 0)
            { // reversing order for being clockwise filling TGeoArb8
                vertices_Cry_Wrap[i][2 * j] = points_cry_local[8 * i + 3 - j].X();
                vertices_Cry_Wrap[i][2 * j + 1] = points_cry_local[8 * i + 3 - j].Y();
            }
            else
            {
                vertices_Cry_Wrap[i][2 * j] = points_cry_local[8 * i + 11 - j].X();
                vertices_Cry_Wrap[i][2 * j + 1] = points_cry_local[8 * i + 11 - j].Y();
            }
        }
    }
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
    {
        for (Int_t j = 0; j < 16; j++)
        {
            if (vertices_Cry_Wrap[i][j] > 0)
                vertices_Cry[i][j] = vertices_Cry_Wrap[i][j] - wrapping_thickness;
            else
                vertices_Cry[i][j] = vertices_Cry_Wrap[i][j] + wrapping_thickness;
        }
    }

    // CEPA USC PART
    //  centers of faces
    TVector3 center_CEPA[N_ALV_TYPES_CEPA * 2]; // 3 geometries, 2 face centers, (3*2)
    TVector3 x_uni_CEPA[N_ALV_TYPES_CEPA * 2];  // unit vectors for each face
    TVector3 y_uni_CEPA[N_ALV_TYPES_CEPA * 2];
    TVector3 z_uni_CEPA[N_ALV_TYPES_CEPA * 2];
    TVector3 center_inn_CEPA[N_ALV_TYPES_CEPA * 2]; // 3 geometries, 2 face centers, (3*2)
    TVector3 x_inn_uni_CEPA[N_ALV_TYPES_CEPA * 2];  // unit vectors for each face
    TVector3 y_inn_uni_CEPA[N_ALV_TYPES_CEPA * 2];
    TVector3 z_inn_uni_CEPA[N_ALV_TYPES_CEPA * 2];
    TVector3 center_cry_CEPA[N_CRY_TYPES_CEPA * 2]; // 12 types of crystals, 2 face centers (12*2)
    TVector3 x_uni_cry_CEPA[N_CRY_TYPES_CEPA * 2];  // unit vectors for each face for crystals
    TVector3 y_uni_cry_CEPA[N_CRY_TYPES_CEPA * 2];
    TVector3 z_uni_cry_CEPA[N_CRY_TYPES_CEPA * 2];
    TRotation rot_CEPA[N_ALV_TYPES_CEPA * 2];     // calculated in each face
    TRotation rot_inn_CEPA[N_ALV_TYPES_CEPA * 2]; // calculated in each face
    TRotation rot_cry_CEPA[N_CRY_TYPES_CEPA * 2]; // only a few are really different if all is ok
    // volume centers
    TVector3 alv_cm_CEPA[N_ALV_TYPES_CEPA];     // 3 geometries
    TVector3 alv_cm_rot_CEPA[N_ALV_TYPES_CEPA]; // 3 geometries, after final rotation
    TVector3 alv_inn_cm_CEPA[N_ALV_TYPES_CEPA]; // 3 geometries
    TVector3 cry_cm_CEPA[N_CRY_TYPES_CEPA];     // 12 types of crystals

    // The center of the crystal faces are first calculated. Then, the unit vectors defining the axis in each faces
    // Third, the rotation moving from the lab system to the unit vectors previously found. To define the
    // volume in Arb8 style, we need the 8 corners in the local frustrum coordinates. Then, we should express
    // the vertices in the coordinate system of the volume center of mass (cm)
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA * 2; i++)
    { // 14 types of crystals, 2 face centers (14*2)
        center_cry_CEPA[i] = points_cry_CEPA[i * 4] + points_cry_CEPA[i * 4 + 1] + points_cry_CEPA[i * 4 + 2] +
                             points_cry_CEPA[i * 4 + 3]; // face centers
        center_cry_CEPA[i] *= 0.25;                      // face centers
        // center_cry_CEPA[i].Print();
        z_uni_cry_CEPA[i] = (points_cry_CEPA[i * 4 + 1] - points_cry_CEPA[i * 4])
                                .Cross(points_cry_CEPA[i * 4 + 2] - points_cry_CEPA[i * 4 + 1]);
        z_uni_cry_CEPA[i] = z_uni_cry_CEPA[i].Unit(); // normal to face center
        x_uni_cry_CEPA[i] =
            points_cry_CEPA[i * 4 + 2] - points_cry_CEPA[i * 4 + 1]; // MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_uni_cry_CEPA[i] = x_uni_cry_CEPA[i].Unit();                   // unit along X
        y_uni_cry_CEPA[i] = z_uni_cry_CEPA[i].Cross(x_uni_cry_CEPA[i]); // unit along Y
        // x_uni_cry_CEPA[i].Print();  y_uni_cry_CEPA[i].Print();    z_uni_cry_CEPA[i].Print();

        // calculate rotation matrix (should be repeated 4 times, just checking)
        rot_cry_CEPA[i].SetZAxis(z_uni_cry_CEPA[i], x_uni_cry_CEPA[i]);
    }
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    { // 14 types of crystals
        cry_cm_CEPA[i] = center_cry_CEPA[2 * i] + center_cry_CEPA[2 * i + 1];
        cry_cm_CEPA[i] *= 0.5; // volume center for each cry
        // cry_cm_CEPA[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each crystal
            points_cry_local_CEPA[i * 8 + j] =
                rot_cry_CEPA[2 * i].Inverse() * (points_cry_CEPA[i * 8 + j] - cry_cm_CEPA[i]);
        }
    }
    /*
        cout<< "Crystal points (local in cm coordinates and global) before reduction: "<< endl;
        for (Int_t i = 0; i < 1; i++)
            for (Int_t j = 0; j < 8; j++){
                points_cry_local_CEPA[i*8+j].Print();
                points_cry_CEPA[i*8+j].Print();
                cout << endl;
            }
        cout<< "END OF Crystal points (local in cm coordinates and global): "<< endl;
    */
    // Reducing the crystals from the UVIGO plans, as they seem to be too large...
    TVector3 reducedis0(crystal_reduction, -crystal_reduction, crystal_reduction);
    TVector3 reducedis1(crystal_reduction, crystal_reduction, crystal_reduction);
    TVector3 reducedis2(-crystal_reduction, crystal_reduction, crystal_reduction);
    TVector3 reducedis3(-crystal_reduction, -crystal_reduction, crystal_reduction);
    TVector3 reducedis4(crystal_reduction, -crystal_reduction, -crystal_reduction);
    TVector3 reducedis5(crystal_reduction, crystal_reduction, -crystal_reduction);
    TVector3 reducedis6(-crystal_reduction, crystal_reduction, -crystal_reduction);
    TVector3 reducedis7(-crystal_reduction, -crystal_reduction, -crystal_reduction);

    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    { // for 14 geometries
        points_cry_local_CEPA[i * 8 + 0] = points_cry_local_CEPA[i * 8 + 0] - reducedis0;
        points_cry_local_CEPA[i * 8 + 1] = points_cry_local_CEPA[i * 8 + 1] - reducedis1;
        points_cry_local_CEPA[i * 8 + 2] = points_cry_local_CEPA[i * 8 + 2] - reducedis2;
        points_cry_local_CEPA[i * 8 + 3] = points_cry_local_CEPA[i * 8 + 3] - reducedis3;
        points_cry_local_CEPA[i * 8 + 4] = points_cry_local_CEPA[i * 8 + 4] - reducedis4;
        points_cry_local_CEPA[i * 8 + 5] = points_cry_local_CEPA[i * 8 + 5] - reducedis5;
        points_cry_local_CEPA[i * 8 + 6] = points_cry_local_CEPA[i * 8 + 6] - reducedis6;
        points_cry_local_CEPA[i * 8 + 7] = points_cry_local_CEPA[i * 8 + 7] - reducedis7;
        // for (Int_t j = 0; j < 8; j++) points_local_CEPA[i*8+j].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices (inverse of the normal point to point local)
            points_cry_CEPA[i * 8 + j] = cry_cm_CEPA[i] + rot_cry_CEPA[2 * i] * points_cry_local_CEPA[i * 8 + j];
            // cout<< "Points in cm coordinates: "<< endl;
            // points_inn_CEPA[i*8+j].Print();
        }
    }
    /*
        cout<< "Crystal points (local in cm coordinates and global) after reduction: "<< endl;
        for (Int_t i = 0; i < 1; i++)
            for (Int_t j = 0; j < 8; j++){
                points_cry_local_CEPA[i*8+j].Print();
                points_cry_CEPA[i*8+j].Print();
                cout << endl;
            }
        cout<< "END OF Crystal points (local in cm coordinates and global): "<< endl;
    */

    // Adding the wrapping
    Double_t partAddition = wrapping_thickness_CEPA;
    TVector3 adddis0(partAddition, -partAddition, partAddition);
    TVector3 adddis1(partAddition, partAddition, partAddition);
    TVector3 adddis2(-partAddition, partAddition, partAddition);
    TVector3 adddis3(-partAddition, -partAddition, partAddition);
    TVector3 adddis4(partAddition, -partAddition, -partAddition);
    TVector3 adddis5(partAddition, partAddition, -partAddition);
    TVector3 adddis6(-partAddition, partAddition, -partAddition);
    TVector3 adddis7(-partAddition, -partAddition, -partAddition);
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    { // for 14 geometries
        points_cry_wrap_local_CEPA[i * 8 + 0] = points_cry_local_CEPA[i * 8 + 0] - adddis0;
        points_cry_wrap_local_CEPA[i * 8 + 1] = points_cry_local_CEPA[i * 8 + 1] - adddis1;
        points_cry_wrap_local_CEPA[i * 8 + 2] = points_cry_local_CEPA[i * 8 + 2] - adddis2;
        points_cry_wrap_local_CEPA[i * 8 + 3] = points_cry_local_CEPA[i * 8 + 3] - adddis3;
        points_cry_wrap_local_CEPA[i * 8 + 4] = points_cry_local_CEPA[i * 8 + 4] - adddis4;
        points_cry_wrap_local_CEPA[i * 8 + 5] = points_cry_local_CEPA[i * 8 + 5] - adddis5;
        points_cry_wrap_local_CEPA[i * 8 + 6] = points_cry_local_CEPA[i * 8 + 6] - adddis6;
        points_cry_wrap_local_CEPA[i * 8 + 7] = points_cry_local_CEPA[i * 8 + 7] - adddis7;
        // for (Int_t j = 0; j < 8; j++) points_local_CEPA[i*8+j].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices (inverse of the normal point to point local)
            points_cry_wrap_CEPA[i * 8 + j] =
                cry_cm_CEPA[i] + rot_cry_CEPA[2 * i] * points_cry_wrap_local_CEPA[i * 8 + j];
            // cout<< "Points in cm coordinates: "<< endl;
            // points_inn_CEPA[i*8+j].Print();
        }
    }
    /*
        cout<< "Crystal_with_wrapping points (local in cm coordinates and global): "<< endl;
        for (Int_t i = 0; i < 1; i++)
            for (Int_t j = 0; j < 8; j++){
                points_cry_wrap_local_CEPA[i*8+j].Print();
                points_cry_wrap_CEPA[i*8+j].Print();
                cout << endl;
            }
        cout<< "END OF Crystal_with_wrapping points (local in cm coordinates and global): "<< endl;
    */

    // internal alveoli corners (calculated crystals corners)
    // They require a later correction to enlarge the alveoli inner space and avoid extrusions
    // as the points_cry_CEPA are given without wrapping
    for (Int_t uOrD = 0; uOrD < 2; uOrD++)
    {                                                                       // upper or lower face
        points_inn_CEPA[uOrD * 4 + 0] = points_cry_wrap_CEPA[uOrD * 4 + 0]; // alveoli 1, first corner of crystal 1
        points_inn_CEPA[uOrD * 4 + 1] = points_cry_wrap_CEPA[1 * 8 + uOrD * 4 + 1]; // second corner of crystal 2
        points_inn_CEPA[uOrD * 4 + 2] = points_cry_wrap_CEPA[2 * 8 + uOrD * 4 + 2]; // third corner of crystal 3
        points_inn_CEPA[uOrD * 4 + 3] = points_cry_wrap_CEPA[3 * 8 + uOrD * 4 + 3]; // fourth corner of crystal 4

        points_inn_CEPA[8 + uOrD * 4 + 0] =
            points_cry_wrap_CEPA[4 * 8 + uOrD * 4 + 0]; // alveoli 2, first corner of crystal 5
        points_inn_CEPA[8 + uOrD * 4 + 1] = points_cry_wrap_CEPA[5 * 8 + uOrD * 4 + 1]; // second corner of crystal 6
        points_inn_CEPA[8 + uOrD * 4 + 2] = points_cry_wrap_CEPA[7 * 8 + uOrD * 4 + 2]; // third corner of crystal 8
        points_inn_CEPA[8 + uOrD * 4 + 3] = points_cry_wrap_CEPA[8 * 8 + uOrD * 4 + 3]; // fourth corner of crystal 9

        points_inn_CEPA[2 * 8 + uOrD * 4 + 0] =
            points_cry_wrap_CEPA[10 * 8 + uOrD * 4 + 0]; // alveoli 3, first corner of crystal 11
        points_inn_CEPA[2 * 8 + uOrD * 4 + 1] =
            points_cry_wrap_CEPA[11 * 8 + uOrD * 4 + 1]; // second corner of crystal 12
        points_inn_CEPA[2 * 8 + uOrD * 4 + 2] =
            points_cry_wrap_CEPA[12 * 8 + uOrD * 4 + 2]; // third corner of crystal 13
        points_inn_CEPA[2 * 8 + uOrD * 4 + 3] =
            points_cry_wrap_CEPA[13 * 8 + uOrD * 4 + 3]; // fourth corner of crystal 14
    }
    /*
        cout<< "Points in cm coordinates: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
            for (Int_t j = 0; j < 8; j++)
                points_inn_CEPA[i*8+j].Print();
        cout<< "END OF Points in cm coordinates: "<< endl;
    */

    // The center of the faces are first calculated. Then, the unit vectors defining the axis in each faces
    // Third, the rotation moving from the lab system to the unit vectors previously found. To define the
    // volume in Arb8 style, we need the 8 corners in the local frustrum coordinates. Then, we should express
    // the vertices in the coordinate system of the volume center of mass (cm)
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA * 2; i++)
    { // for 3 geometries, 2 face centers (3*2)
        center_inn_CEPA[i] = points_inn_CEPA[i * 4] + points_inn_CEPA[i * 4 + 1] + points_inn_CEPA[i * 4 + 2] +
                             points_inn_CEPA[i * 4 + 3]; // face centers
        center_inn_CEPA[i] *= 0.25;                      // face centers
        // cout<< "Alv centers coordinates: "<< endl; center_inn_CEPA[i].Print(); cout<< "END OF Alv centers
        // coordinates: "<< endl;
        z_inn_uni_CEPA[i] = (points_inn_CEPA[i * 4 + 1] - points_inn_CEPA[i * 4])
                                .Cross(points_inn_CEPA[i * 4 + 2] - points_inn_CEPA[i * 4 + 1]);
        z_inn_uni_CEPA[i] = z_inn_uni_CEPA[i].Unit(); // normal to face center
        x_inn_uni_CEPA[i] =
            points_inn_CEPA[i * 4 + 2] - points_inn_CEPA[i * 4 + 1]; // MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_inn_uni_CEPA[i] = x_inn_uni_CEPA[i].Unit();                   // unit along X
        y_inn_uni_CEPA[i] = z_inn_uni_CEPA[i].Cross(x_inn_uni_CEPA[i]); // unit along Y
        // x_inn_uni_CEPA[i].Print();  y_inn_uni_CEPA[i].Print();   z_inn_uni_CEPA[i].Print();
        //  calculate rotation matrix for the 3 geometries (should be repeated 4 times, just checking)
        rot_inn_CEPA[i].SetZAxis(z_inn_uni_CEPA[i], x_inn_uni_CEPA[i]);
    }
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    { // for 3 geometries
        alv_inn_cm_CEPA[i] = center_inn_CEPA[2 * i] + center_inn_CEPA[2 * i + 1];
        alv_inn_cm_CEPA[i] *= 0.5; // volume center for all alv (outer and inner)
        // alv_inn_cm_CEPA[i].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each alveolus
            points_inn_local_CEPA[i * 8 + j] =
                rot_inn_CEPA[2 * i].Inverse() * (points_inn_CEPA[i * 8 + j] - alv_inn_cm_CEPA[i]);
            // cout<< "Points in cm coordinates: "<< endl;
            // points_inn_local_CEPA[i*8+j].Print();
        }
    }

    // Adding now the wrapping and security margin to the inner alveoli points
    Double_t totalAddition = security_margin;
    TVector3 dis0(totalAddition, -totalAddition, totalAddition);
    TVector3 dis1(totalAddition, totalAddition, totalAddition);
    TVector3 dis2(-totalAddition, totalAddition, totalAddition);
    TVector3 dis3(-totalAddition, -totalAddition, totalAddition);
    TVector3 dis4(totalAddition, -totalAddition, -totalAddition);
    TVector3 dis5(totalAddition, totalAddition, -totalAddition);
    TVector3 dis6(-totalAddition, totalAddition, -totalAddition);
    TVector3 dis7(-totalAddition, -totalAddition, -totalAddition);
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    { // for 3 geometries
        points_inn_local_CEPA[i * 8 + 0] = points_inn_local_CEPA[i * 8 + 0] - dis0;
        points_inn_local_CEPA[i * 8 + 1] = points_inn_local_CEPA[i * 8 + 1] - dis1;
        points_inn_local_CEPA[i * 8 + 2] = points_inn_local_CEPA[i * 8 + 2] - dis2;
        points_inn_local_CEPA[i * 8 + 3] = points_inn_local_CEPA[i * 8 + 3] - dis3;
        points_inn_local_CEPA[i * 8 + 4] = points_inn_local_CEPA[i * 8 + 4] - dis4;
        points_inn_local_CEPA[i * 8 + 5] = points_inn_local_CEPA[i * 8 + 5] - dis5;
        points_inn_local_CEPA[i * 8 + 6] = points_inn_local_CEPA[i * 8 + 6] - dis6;
        points_inn_local_CEPA[i * 8 + 7] = points_inn_local_CEPA[i * 8 + 7] - dis7;
        // for (Int_t j = 0; j < 8; j++) points_local_CEPA[i*8+j].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each inner alveolus (inverse of the normal point to point local)
            points_inn_CEPA[i * 8 + j] = alv_inn_cm_CEPA[i] + rot_inn_CEPA[2 * i] * points_inn_local_CEPA[i * 8 + j];
            // cout<< "Points in cm coordinates: "<< endl;
            // points_inn_CEPA[i*8+j].Print();
        }
    }
    /*
        cout<< "Points (inner) in cm coordinates: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
            for (Int_t j = 0; j < 8; j++){
                points_inn_local_CEPA[i*8+j].Print();
                points_inn_CEPA[i*8+j].Print();
            }
    */

    // calculating the OUTER points from the INNER ones
    TVector3 d0(cf_thickness_CEPA, -cf_thickness_CEPA, cf_thickness_CEPA);
    TVector3 d1(cf_thickness_CEPA, cf_thickness_CEPA, cf_thickness_CEPA);
    TVector3 d2(-cf_thickness_CEPA, cf_thickness_CEPA, cf_thickness_CEPA);
    TVector3 d3(-cf_thickness_CEPA, -cf_thickness_CEPA, cf_thickness_CEPA);
    TVector3 d4(cf_thickness_CEPA, -cf_thickness_CEPA, -cf_thickness_CEPA);
    TVector3 d5(cf_thickness_CEPA, cf_thickness_CEPA, -cf_thickness_CEPA);
    TVector3 d6(-cf_thickness_CEPA, cf_thickness_CEPA, -cf_thickness_CEPA);
    TVector3 d7(-cf_thickness_CEPA, -cf_thickness_CEPA, -cf_thickness_CEPA);
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    { // for 3 geometries
        points_local_CEPA[i * 8 + 0] = points_inn_local_CEPA[i * 8 + 0] - d0;
        points_local_CEPA[i * 8 + 1] = points_inn_local_CEPA[i * 8 + 1] - d1;
        points_local_CEPA[i * 8 + 2] = points_inn_local_CEPA[i * 8 + 2] - d2;
        points_local_CEPA[i * 8 + 3] = points_inn_local_CEPA[i * 8 + 3] - d3;
        points_local_CEPA[i * 8 + 4] = points_inn_local_CEPA[i * 8 + 4] - d4;
        points_local_CEPA[i * 8 + 5] = points_inn_local_CEPA[i * 8 + 5] - d5;
        points_local_CEPA[i * 8 + 6] = points_inn_local_CEPA[i * 8 + 6] - d6;
        points_local_CEPA[i * 8 + 7] = points_inn_local_CEPA[i * 8 + 7] - d7;
        // for (Int_t j = 0; j < 8; j++) points_local_CEPA[i*8+j].Print();
        for (Int_t j = 0; j < 8; j++)
        { // for the 8 vertices of each inner alveolus (inverse of the normal point to point local)
            points_CEPA[i * 8 + j] = alv_inn_cm_CEPA[i] + rot_inn_CEPA[2 * i] * points_local_CEPA[i * 8 + j];
            // cout<< "Points in cm coordinates: "<< endl;
            // points_inn_CEPA[i*8+j].Print();
        }
    }
    /*
        cout<< "Points (outer) in cm coordinates: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
            for (Int_t j = 0; j < 8; j++){
                points_local_CEPA[i*8+j].Print();
                points_CEPA[i*8+j].Print();
            }
    */

    for (Int_t i = 0; i < N_ALV_TYPES_CEPA * 2; i++)
    { // for 3 geometries, 2 face centers (3*2)
        center_CEPA[i] = points_CEPA[i * 4] + points_CEPA[i * 4 + 1] + points_CEPA[i * 4 + 2] +
                         points_CEPA[i * 4 + 3]; // face centers
        center_CEPA[i] *= 0.25;                  // face centers
        // center_CEPA[i].Print();
        z_uni_CEPA[i] =
            (points_CEPA[i * 4 + 1] - points_CEPA[i * 4]).Cross(points_CEPA[i * 4 + 2] - points_CEPA[i * 4 + 1]);
        z_uni_CEPA[i] = z_uni_CEPA[i].Unit();                            // normal to face center
        x_uni_CEPA[i] = points_CEPA[i * 4 + 2] - points_CEPA[i * 4 + 1]; // MODIFIED FROM BARREL+IPHOS DEFINITION!!!!

        x_uni_CEPA[i] = x_uni_CEPA[i].Unit();               // unit along X
        y_uni_CEPA[i] = z_uni_CEPA[i].Cross(x_uni_CEPA[i]); // unit along Y
        // x_uni_CEPA[i].Print();  y_uni_CEPA[i].Print();   z_uni_CEPA[i].Print();
        // calculate rotation matrix for the 3 geometries (should be repeated 4 times, just checking)
        rot_CEPA[i].SetZAxis(z_uni_CEPA[i], x_uni_CEPA[i]);
    }

    // outer alv
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    { // for 3 geometries
        alv_cm_CEPA[i] = center_CEPA[2 * i] + center_CEPA[2 * i + 1];
        alv_cm_CEPA[i] *= 0.5; // volume center for all alv (outer and inner)
        // alv_cm_CEPA[i].Print();
        // NEXT LINES ARE ONLY FOR CHECKING THE ROTATION
        // for (Int_t j = 0; j < 8; j++)
        //{ // for the 8 vertices of each alveolus
        // points_CEPA_local_CEPA[i * 8 + j] = rot_CEPA[2 * i].Inverse() * (points_CEPA[i * 8 + j] - alv_cm_CEPA[i]);
        // cout<< "Points in cm coordinates: "<< endl;
        // points_local_CEPA[i*8+j].Print();
        //}
    }

    // location of the crystals in the alveoli
    TVector3 cry_position_CEPA[N_CRY_TYPES_CEPA];       // 14 types of crystals
    TVector3 cry_position_local_CEPA[N_CRY_TYPES_CEPA]; // 14 types of crystals

    // Relative Crystal rotation in each alveoli. Obtained from the crystal unit vector in
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    { // four crystals per alv
        if (i < 4)
        { // first alveoli crystals
            cry_position_CEPA[i] = cry_cm_CEPA[i] - alv_cm_CEPA[0];
            cry_position_local_CEPA[i] = rot_CEPA[0].Inverse() * (cry_cm_CEPA[i] - alv_cm_CEPA[0]);
        }
        if (i > 3 && i < 10)
        { // second alveoli crystals
            cry_position_CEPA[i] = cry_cm_CEPA[i] - alv_cm_CEPA[1];
            cry_position_local_CEPA[i] = rot_CEPA[2].Inverse() * (cry_cm_CEPA[i] - alv_cm_CEPA[1]);
        }
        if (i > 9)
        {
            cry_position_CEPA[i] = cry_cm_CEPA[i] - alv_cm_CEPA[2]; // third alveoli crystals
            cry_position_local_CEPA[i] = rot_CEPA[4].Inverse() * (cry_cm_CEPA[i] - alv_cm_CEPA[2]);
        }
    }
    /*
            cout<< "Points local in cm coordinates: "<< endl;
            for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
                for (Int_t j = 0; j < 8; j++)
                     points_local_CEPA[i*8+j].Print();
            cout<< "Points inn local in cm coordinates: "<< endl;
            for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
                for (Int_t j = 0; j < 8; j++)
                    points_inn_local_CEPA[i*8+j].Print();
            cout<< "Points (check) cry local in cm coordinates after reduction: "<< endl;
            for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
                for (Int_t j = 0; j < 8; j++)
                    points_cry_local_CEPA[i * 8 + j].Print();
    */

    // Redefinition of vertices for the construction of the Alveoli, using TGeoArb8
    Double_t* vertices_Alv_CEPA[N_ALV_TYPES_CEPA]; // 3 geometries
    Double_t* vertices_inner_Alv_CEPA[N_ALV_TYPES_CEPA];
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    {
        vertices_Alv_CEPA[i] = new Double_t[16];
        vertices_inner_Alv_CEPA[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        { // reversing order for being clockwise filling TGeoArb8
            if (j == 1 || j == 5)
            {
                vertices_Alv_CEPA[i][2 * j] = points_local_CEPA[8 * i + j + 2].X();
                vertices_Alv_CEPA[i][2 * j + 1] = points_local_CEPA[8 * i + j + 2].Y();
                vertices_inner_Alv_CEPA[i][2 * j] = points_inn_local_CEPA[8 * i + j + 2].X();
                vertices_inner_Alv_CEPA[i][2 * j + 1] = points_inn_local_CEPA[8 * i + j + 2].Y();
            }
            else if (j == 3 || j == 7)
            {
                vertices_Alv_CEPA[i][2 * j] = points_local_CEPA[8 * i + j - 2].X();
                vertices_Alv_CEPA[i][2 * j + 1] = points_local_CEPA[8 * i + j - 2].Y();
                vertices_inner_Alv_CEPA[i][2 * j] = points_inn_local_CEPA[8 * i + j - 2].X();
                vertices_inner_Alv_CEPA[i][2 * j + 1] = points_inn_local_CEPA[8 * i + j - 2].Y();
            }
            else
            {
                vertices_Alv_CEPA[i][2 * j] = points_local_CEPA[8 * i + j].X();
                vertices_Alv_CEPA[i][2 * j + 1] = points_local_CEPA[8 * i + j].Y();
                vertices_inner_Alv_CEPA[i][2 * j] = points_inn_local_CEPA[8 * i + j].X();
                vertices_inner_Alv_CEPA[i][2 * j + 1] = points_inn_local_CEPA[8 * i + j].Y();
            }
        }
    }
    /*
        cout<< "Points local in cm coordinates reordered: "<< endl;
        for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
        {
            for (Int_t j = 0; j < 16; j++)
            {
            cout << vertices_Alv_CEPA[i][j] << endl;
            }
        }
    */

    // Redefinition of vertices for the construction of the Crystals, using TGeoArb8
    // In this version, points correspond to the crystal and the wrapping should be added
    Double_t* vertices_Cry_CEPA[N_CRY_TYPES_CEPA];      // (14 crystal types)
    Double_t* vertices_Cry_Wrap_CEPA[N_CRY_TYPES_CEPA]; // (14 crystal types)
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    {
        vertices_Cry_CEPA[i] = new Double_t[16];
        vertices_Cry_Wrap_CEPA[i] = new Double_t[16];
    }
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
            if ((4 - j) > 0)
            { // reversing order for being clockwise filling TGeoArb8
                vertices_Cry_CEPA[i][2 * j] = points_cry_local_CEPA[8 * i + 3 - j].X();
                vertices_Cry_CEPA[i][2 * j + 1] = points_cry_local_CEPA[8 * i + 3 - j].Y();
                vertices_Cry_Wrap_CEPA[i][2 * j] = points_cry_wrap_local_CEPA[8 * i + 3 - j].X();
                vertices_Cry_Wrap_CEPA[i][2 * j + 1] = points_cry_wrap_local_CEPA[8 * i + 3 - j].Y();
            }
            else
            {
                vertices_Cry_CEPA[i][2 * j] = points_cry_local_CEPA[8 * i + 11 - j].X();
                vertices_Cry_CEPA[i][2 * j + 1] = points_cry_local_CEPA[8 * i + 11 - j].Y();
                vertices_Cry_Wrap_CEPA[i][2 * j] = points_cry_wrap_local_CEPA[8 * i + 11 - j].X();
                vertices_Cry_Wrap_CEPA[i][2 * j + 1] = points_cry_wrap_local_CEPA[8 * i + 11 - j].Y();
            }
        }
    }

    // CONSTRUCTING BARREL+iPHOS PART
    TGeoVolume** Alv_vol;
    Alv_vol = new TGeoVolume*[N_ALV_TYPES];
    TGeoVolume** Alv_inner_vol;
    Alv_inner_vol = new TGeoVolume*[N_ALV_TYPES];
    TGeoVolume** Cry_vol_wrap;
    Cry_vol_wrap = new TGeoVolume*[N_CRY_TYPES];
    TGeoVolume** Cry_vol;
    Cry_vol = new TGeoVolume*[N_CRY_TYPES];

    TString AlvGlobalName = "Alveolus_";
    TString AlvGlobalNameInner = "InnerAlv_";
    // TString name_Alv[N_ALV_TYPES] = { "BR_16a", "BR_15a", "BR_14a", "BR_13a", "BR_12a", "BR_11a", "BR_10a", "BR_09a",
    //                                  "BR_08a", "BR_07a", "BR_06a", "BR_05a", "BR_04a", "BR_03a", "BR_02a", "BR_01a",
    //                                  "iP_04a", "iP_03a", "iP_02a", "iP_01d", "iP_01c", "iP_01b", "iP_01a" };
    // Substitute names in previous array (CAD names) to simplify the R3BRoot code
    TString name_Alv[N_ALV_TYPES] = { "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12",
                                      "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23" };
    Double_t halfLengthAlv[N_ALV_TYPES] = { 8.5,  9.5,  9.5,  9.5,  9.5,  10.5, 10.5, 11.0, 11.0, 11.0, 11.5, 11.5,
                                            11.5, 13.5, 13.5, 13.5, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6, 13.6 }; // cm
    Double_t halfLengthAlv_inner[N_ALV_TYPES] = { 8.485,  9.485,  9.485,  9.485,  9.485,  10.485, 10.485, 10.985,
                                                  10.985, 10.985, 11.485, 11.485, 11.485, 13.485, 13.485, 13.485,
                                                  13.585, 13.585, 13.585, 13.585, 13.585, 13.585, 13.585 }; // cm
    TString WrapCryGlobalName = "WrapCry_";
    TString CryGlobalName = "Crystal_";
    TString name_Cry[4] = { "_1", "_2", "_3", "_4" };
    // TODO, KNOWN ISSUE TO SOLVE!
    // Crystal BR_16_1 has the same length and size as the inner alv (-1mm) (no crystal was defined in CAD)
    Double_t halfLengthCry[N_ALV_TYPES] = { 8.385,   7.0,     7.0,     7.0,    7.0,     8.0,     8.0,   8.5,
                                            8.5,     8.5,     9.0,     9.0,    9.0,     11.0,    11.0,  11.0,
                                            10.9754, 10.9799, 10.9841, 10.985, 10.9818, 10.9819, 10.983 }; // cm

    // length checks
    /*
    cout << "Check alveoli length (below 1 microns). No error message indicates that all is ok..." << endl;
    for(Int_t i=0;i<N_ALV_TYPES;i++) cout << halfLengthAlv[i] << "\t";
    cout << endl;
    for(Int_t i=0;i<N_ALV_TYPES;i++) cout << halfLengthAlv_inner[i] << "\t";
    cout << endl;
    for(Int_t i=0;i<N_ALV_TYPES;i++){
      for(Int_t j=0;j<8;j++){
        if(abs(abs(points_local[i*16+j].Z())-halfLengthAlv[i]) > 0.0001) {
      cout << "error! Check alv problem at i=" << i << ", j=" << j<< endl;
      cout << abs(points_local[i*8+j].Z()) << " " << halfLengthAlv[i] << endl;
        }
        if(abs(abs(points_local[i*16+j+8].Z())-halfLengthAlv_inner[i]) > 0.0001) {
      cout << "error! Check inner alv problem at i=" << i << ", j=" << j<< endl;
      cout << abs(points_local[i*8+j+8].Z()) << " " << halfLengthAlv_inner[i] << endl;
        }
      }
    }
    cout << "Check crystal length (deviation below 1 micron). No error message indicates that all is ok..." << endl;
    for(Int_t i=0;i<N_ALV_TYPES;i++) cout << halfLengthCry[i] << "\t";
    cout << endl;
    for(Int_t i=0;i<N_ALV_TYPES;i++){
      for(Int_t k=0;k<4;k++){ //four crystal per alveoli
        for(Int_t j=0;j<8;j++){
      //cout << "  i=" << i << ",  k=" << k <<  ",  j=" << j <<  ",  a=" << i*8+j
      //     <<",  b=" << (i-1)*32+k*8+j+8 << ",  c=" << (i-19)*24+k*8+j+584 << endl;
      if(i== 0 && k==0 && abs(abs(points_cry_local[i*8+j].Z())-halfLengthCry[i]) > 0.0001) {
        cout << "error! Check crystal problem at i=0" << ", j=" << j << endl;
        cout << abs(points_cry_local[i*8+j].Z()) << " " << halfLengthCry[i] << endl;
      }
      else if(i> 0 && i< 19 && abs(abs(points_cry_local[(i-1)*32+k*8+j+8].Z())-halfLengthCry[i]) > 0.0001) {
        cout << "error! Check crystal problem at i=" << i << ", j=" << j << ", k=" << k << endl;
        cout << abs(points_cry_local[(i-1)*32+k*8+j+8].Z()) << " " << halfLengthCry[i] << endl;
      }
      else if(i>18 && k<3 && abs(abs(points_cry_local[(i-19)*24+k*8+j+584].Z())-halfLengthCry[i]) > 0.0001) {
        cout << "error! Check crystal problem at i=" << i << ", j=" << j << ", k=" << k << endl;
        cout << abs(points_cry_local[(i-19)*24+k*8+j+584].Z()) << " " << halfLengthCry[i] << endl;
      }
        }
      }
    }
    */

    TGeoRotation* rotUni = new TGeoRotation();
    TGeoRotation** rotAlv = new TGeoRotation*[N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
        rotAlv[i] = new TGeoRotation();
    TGeoRotation** rotCry = new TGeoRotation*[N_CRY_TYPES];
    for (Int_t i = 0; i < N_CRY_TYPES; i++)
        rotCry[i] = new TGeoRotation();
    Double_t rotEle[9];

    // rotation
    TGeoRotation** rotOnZ = new TGeoRotation*[32];
    for (Int_t i = 0; i < 32; i++)
        rotOnZ[i] = new TGeoRotation();
    for (Int_t i = 0; i < 32; i++)
        rotOnZ[i]->RotateZ(-11.25 * i);
    TRotation** rotationOnZ = new TRotation*[32];
    for (Int_t i = 0; i < 32; i++)
        rotationOnZ[i] = new TRotation();
    for (Int_t i = 0; i < 32; i++)
        rotationOnZ[i]->RotateZ(i * (-11.25 * TMath::Pi()) / 180);

    TGeoRotation** rotAlvFinal = new TGeoRotation*[32 * N_ALV_TYPES];
    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        for (Int_t j = 0; j < 32; j++)
        {
            rotAlvFinal[i * 32 + j] = new TGeoRotation((*rotOnZ[j]) * (*rotAlv[i]));
        }
    }

    for (Int_t i = 0; i < N_ALV_TYPES; i++)
    {
        Alv_vol[i] =
            gGeoManager->MakeArb8(AlvGlobalName + name_Alv[i], pCarbonFibreMedium, halfLengthAlv[i], vertices_Alv[i]);
        Alv_vol[i]->SetLineColor(kBlue);
        Alv_vol[i]->SetVisLeaves(kTRUE);
        Alv_vol[i]->SetVisibility(kTRUE);
        Alv_vol[i]->SetVisContainers(kTRUE);

        Alv_inner_vol[i] = gGeoManager->MakeArb8(
            AlvGlobalNameInner + name_Alv[i], pAirMedium, halfLengthAlv_inner[i], vertices_inner_Alv[i]);
        Alv_inner_vol[i]->SetLineColor(kRed);
        Alv_inner_vol[i]->SetVisLeaves(kTRUE);
        Alv_inner_vol[i]->SetVisibility(kTRUE);
        Alv_inner_vol[i]->SetVisContainers(kTRUE);

        if (i == 0)
        { // first alveolus with a single crystal
            Cry_vol[i] = gGeoManager->MakeArb8(CryGlobalName + name_Alv[i] + name_Cry[i],
                                               pCsIMedium,
                                               halfLengthCry[i] - wrapping_thickness,
                                               vertices_Cry[i]);
            Cry_vol[i]->SetLineColor(kMagenta);
            Cry_vol[i]->SetVisLeaves(kTRUE);
            Cry_vol[i]->SetVisibility(kTRUE);
            Cry_vol[i]->SetVisContainers(kTRUE);

            Cry_vol_wrap[i] = gGeoManager->MakeArb8(
                WrapCryGlobalName + name_Alv[i] + name_Cry[i], pWrappingMedium, halfLengthCry[i], vertices_Cry_Wrap[i]);
            Cry_vol_wrap[i]->SetLineColor(kGreen);
            Cry_vol_wrap[i]->SetVisLeaves(kTRUE);
            Cry_vol_wrap[i]->SetVisibility(kTRUE);
            Cry_vol_wrap[i]->SetVisContainers(kTRUE);

            Cry_vol_wrap[i]->AddNode(Cry_vol[i], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
            Alv_inner_vol[i]->AddNode(
                Cry_vol_wrap[i],
                0,
                new TGeoCombiTrans(
                    cry_position_local[i].X(), cry_position_local[i].Y(), cry_position_local[i].Z(), rotUni));
        }
        else if (i < 19)
        { // four crystals per alv
            for (Int_t j = 0; j < 4; j++)
            {
                Cry_vol[4 * i - 3 + j] = gGeoManager->MakeArb8(CryGlobalName + name_Alv[i] + name_Cry[j],
                                                               pCsIMedium,
                                                               halfLengthCry[i] - wrapping_thickness,
                                                               vertices_Cry[4 * i - 3 + j]);
                Cry_vol[4 * i - 3 + j]->SetLineColor(kMagenta);
                Cry_vol[4 * i - 3 + j]->SetVisLeaves(kTRUE);
                Cry_vol[4 * i - 3 + j]->SetVisibility(kTRUE);
                Cry_vol[4 * i - 3 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap[4 * i - 3 + j] = gGeoManager->MakeArb8(WrapCryGlobalName + name_Alv[i] + name_Cry[j],
                                                                    pWrappingMedium,
                                                                    halfLengthCry[i],
                                                                    vertices_Cry_Wrap[4 * i - 3 + j]);
                Cry_vol_wrap[4 * i - 3 + j]->SetLineColor(kGreen);
                Cry_vol_wrap[4 * i - 3 + j]->SetVisLeaves(kTRUE);
                Cry_vol_wrap[4 * i - 3 + j]->SetVisibility(kTRUE);
                Cry_vol_wrap[4 * i - 3 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap[4 * i - 3 + j]->AddNode(Cry_vol[4 * i - 3 + j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
                Alv_inner_vol[i]->AddNode(Cry_vol_wrap[4 * i - 3 + j],
                                          0,
                                          new TGeoCombiTrans(cry_position_local[4 * i - 3 + j].X(),
                                                             cry_position_local[4 * i - 3 + j].Y(),
                                                             cry_position_local[4 * i - 3 + j].Z(),
                                                             rotUni));
            }
        }
        else
        { // last four alveoli with three crystals per alv
            for (Int_t j = 0; j < 3; j++)
            {
                Cry_vol[73 + 3 * (i - 19) + j] = gGeoManager->MakeArb8(CryGlobalName + name_Alv[i] + name_Cry[j],
                                                                       pCsIMedium,
                                                                       halfLengthCry[i] - wrapping_thickness,
                                                                       vertices_Cry[73 + 3 * (i - 19) + j]);
                Cry_vol[73 + 3 * (i - 19) + j]->SetLineColor(kMagenta);
                Cry_vol[73 + 3 * (i - 19) + j]->SetVisLeaves(kTRUE);
                Cry_vol[73 + 3 * (i - 19) + j]->SetVisibility(kTRUE);
                Cry_vol[73 + 3 * (i - 19) + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap[73 + 3 * (i - 19) + j] =
                    gGeoManager->MakeArb8(WrapCryGlobalName + name_Alv[i] + name_Cry[j],
                                          pWrappingMedium,
                                          halfLengthCry[i],
                                          vertices_Cry_Wrap[73 + 3 * (i - 19) + j]);
                Cry_vol_wrap[73 + 3 * (i - 19) + j]->SetLineColor(kGreen);
                Cry_vol_wrap[73 + 3 * (i - 19) + j]->SetVisLeaves(kTRUE);
                Cry_vol_wrap[73 + 3 * (i - 19) + j]->SetVisibility(kTRUE);
                Cry_vol_wrap[73 + 3 * (i - 19) + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap[73 + 3 * (i - 19) + j]->AddNode(
                    Cry_vol[73 + 3 * (i - 19) + j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
                Alv_inner_vol[i]->AddNode(Cry_vol_wrap[73 + 3 * (i - 19) + j],
                                          0,
                                          new TGeoCombiTrans(cry_position_local[73 + 3 * (i - 19) + j].X(),
                                                             cry_position_local[73 + 3 * (i - 19) + j].Y(),
                                                             cry_position_local[73 + 3 * (i - 19) + j].Z(),
                                                             rotUni));
            }
        }
        // Inner volume center is displaced 150 microns along Z
        Alv_vol[i]->AddNode(Alv_inner_vol[i], 0, new TGeoCombiTrans(0, 0, 0.015, rotUni));

        rotEle[0] = rot[4 * i].XX();
        rotEle[1] = rot[4 * i].XY();
        rotEle[2] = rot[4 * i].XZ();
        rotEle[3] = rot[4 * i].YX();
        rotEle[4] = rot[4 * i].YY();
        rotEle[5] = rot[4 * i].YZ();
        rotEle[6] = rot[4 * i].ZX();
        rotEle[7] = rot[4 * i].ZY();
        rotEle[8] = rot[4 * i].ZZ();
        rotAlv[i]->SetMatrix(rotEle);

        for (Int_t j = 0; j < 32; j++)
        { // rotation around Z
            // rotAlvFinal[i * 32 + j] = new TGeoRotation((*rotOnZ[j]) * (*rotAlv[i]));
            // alv_cm_rot[2 * i] = (*rotationOnZ[j]) * alv_cm[2 * i];

            // ALIGNMENT DISPLACEMENT FOR SOME EXPERIMENTS
            Double_t disp_halfBarrel = 0;
            if (strncmp(expNumber, "s455", 4) == 0 || strncmp(expNumber, "s515", 4) == 0)
            {
                if (j < 16)
                    disp_halfBarrel = 1.25;
                else
                    disp_halfBarrel = -1.25;
            }

            if (i > 18 && j > 7)
                continue;

            if (i > 18)
            {
                // ALIGNMENT DISPLACEMENT FOR SOME EXPERIMENTS
                if (strncmp(expNumber, "s455", 4) == 0 || strncmp(expNumber, "s515", 4) == 0)
                {
                    if (j < 4)
                        disp_halfBarrel = 1.25;
                    else
                        disp_halfBarrel = -1.25;
                }

                rotAlvFinal[i * 32 + 4 * j] = new TGeoRotation((*rotOnZ[4 * j]) * (*rotAlv[i]));
                alv_cm_rot[2 * i] = (*rotationOnZ[4 * j]) * alv_cm[2 * i];
                if (isCrystalInstalled(i + 1, j, installedCrystals))
                {
                    pWorld->AddNode(Alv_vol[i],
                                    j,
                                    new TGeoCombiTrans(alv_cm_rot[2 * i].X() + disp_halfBarrel,
                                                       alv_cm_rot[2 * i].Y(),
                                                       alv_cm_rot[2 * i].Z(),
                                                       rotAlvFinal[i * 32 + 4 * j]));
                }
            }
            else
            {
                rotAlvFinal[i * 32 + j] = new TGeoRotation((*rotOnZ[j]) * (*rotAlv[i]));
                alv_cm_rot[2 * i] = (*rotationOnZ[j]) * alv_cm[2 * i];
                if (isCrystalInstalled(i + 1, j, installedCrystals))
                {
                    pWorld->AddNode(Alv_vol[i],
                                    j,
                                    new TGeoCombiTrans(alv_cm_rot[2 * i].X() + disp_halfBarrel,
                                                       alv_cm_rot[2 * i].Y(),
                                                       alv_cm_rot[2 * i].Z(),
                                                       rotAlvFinal[i * 32 + j]));
                }
            }
        }
    }

    // CONSTRUCTING CEPA_USC PART
    TGeoVolume** Alv_vol_CEPA;
    Alv_vol_CEPA = new TGeoVolume*[N_ALV_TYPES_CEPA];
    TGeoVolume** Alv_inner_vol_CEPA;
    Alv_inner_vol_CEPA = new TGeoVolume*[N_ALV_TYPES_CEPA];
    TGeoVolume** Cry_vol_wrap_CEPA;
    Cry_vol_wrap_CEPA = new TGeoVolume*[N_CRY_TYPES_CEPA];
    TGeoVolume** Cry_vol_CEPA;
    Cry_vol_CEPA = new TGeoVolume*[N_CRY_TYPES_CEPA];

    TString AlvGlobalName_CEPA = "Alveolus_CCSI_";
    TString AlvGlobalNameInner_CEPA = "InnerAlv_CCSI_";
    // Substitute names in previous array (CAD names) to simplify the R3BRoot code
    TString name_Alv_CEPA[N_ALV_TYPES_CEPA] = { "01", "02", "03" };
    Double_t halfLengthAlv_CEPA[N_ALV_TYPES] = { points_local_CEPA[4].Z(),
                                                 points_local_CEPA[12].Z(),
                                                 points_local_CEPA[20].Z() }; // cm
    Double_t halfLengthAlv_inner_CEPA[N_ALV_TYPES] = { points_inn_local_CEPA[4].Z(),
                                                       points_inn_local_CEPA[12].Z(),
                                                       points_inn_local_CEPA[20].Z() }; // cm

    TString WrapCryGlobalName_CEPA = "WrapCry_CCSI_";
    TString CryGlobalName_CEPA = "Crystal_CCSI_";
    TString name_Cry_CEPA[4] = { "_1", "_2", "_3", "_4" };
    TString name_Cry_CEPA_6[6] = { "_1", "_2", "_3", "_4", "_5", "_6" };
    // For the moment same length as the inner alveoli
    Double_t halfLengthCry_CEPA[N_ALV_TYPES_CEPA] = { points_cry_local_CEPA[4].Z(),
                                                      points_cry_local_CEPA[12].Z(),
                                                      points_cry_local_CEPA[20].Z() }; // cm

    TGeoRotation** rotAlv_CEPA = new TGeoRotation*[N_ALV_TYPES_CEPA];
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
        rotAlv_CEPA[i] = new TGeoRotation();
    TGeoRotation** rotCry_CEPA = new TGeoRotation*[N_CRY_TYPES_CEPA];
    for (Int_t i = 0; i < N_CRY_TYPES_CEPA; i++)
        rotCry_CEPA[i] = new TGeoRotation();
    Double_t rotEle_CEPA[9];

    // rotation
    TGeoRotation** rotOnZ_CEPA = new TGeoRotation*[8];
    for (Int_t i = 0; i < 8; i++)
    {
        rotOnZ_CEPA[i] = new TGeoRotation();
    }
    for (Int_t i = 0; i < 8; i++)
    {
        rotOnZ_CEPA[i]->RotateZ(-45.0 * i); // 67.5 is the offset to put the first alveoli below the first of the barrel
        // rotOnZ[i]->Print();
    }
    TRotation** rotationOnZ_CEPA = new TRotation*[8];
    for (Int_t i = 0; i < 8; i++)
    {
        rotationOnZ_CEPA[i] = new TRotation();
    }
    for (Int_t i = 0; i < 8; i++)
    {
        rotationOnZ_CEPA[i]->RotateZ((i * -45.0) * TMath::Pi() / 180);
    }

    TGeoRotation** rotAlvFinal_CEPA = new TGeoRotation*[8 * N_ALV_TYPES_CEPA];
    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    {
        for (Int_t j = 0; j < 8; j++)
        {
            // rotAlvFinal_CEPA[i * 8 + j] = new TGeoRotation( (*rotOnZ_CEPA[j]) * (*rotAlv_CEPA[i]));
            rotAlvFinal_CEPA[i * 8 + j] = new TGeoRotation((*rotOnZ_CEPA[j]) * (*rotAlv_CEPA[i]));
            // if(i==0) rotAlvFinal_CEPA[i * 8 + j]->Print();
        }
    }

    for (Int_t i = 0; i < N_ALV_TYPES_CEPA; i++)
    {
        Alv_vol_CEPA[i] = gGeoManager->MakeArb8(
            AlvGlobalName_CEPA + name_Alv_CEPA[i], pCarbonFibreMedium, halfLengthAlv_CEPA[i], vertices_Alv_CEPA[i]);
        Alv_vol_CEPA[i]->SetLineColor(kBlue);
        Alv_vol_CEPA[i]->SetVisLeaves(kTRUE);
        Alv_vol_CEPA[i]->SetVisibility(kTRUE);
        Alv_vol_CEPA[i]->SetVisContainers(kTRUE);

        Alv_inner_vol_CEPA[i] = gGeoManager->MakeArb8(AlvGlobalNameInner_CEPA + name_Alv_CEPA[i],
                                                      pAirMedium,
                                                      halfLengthAlv_inner_CEPA[i],
                                                      vertices_inner_Alv_CEPA[i]);
        Alv_inner_vol_CEPA[i]->SetLineColor(kRed);
        Alv_inner_vol_CEPA[i]->SetVisLeaves(kTRUE);
        Alv_inner_vol_CEPA[i]->SetVisibility(kTRUE);
        Alv_inner_vol_CEPA[i]->SetVisContainers(kTRUE);

        // four crystals per alv, but alveolus 2 with 6 crystals
        if (i == 0) // 4 crystal alveoli
        {
            for (Int_t j = 0; j < 4; j++)
            {
                Cry_vol_CEPA[j] = gGeoManager->MakeArb8(CryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA[j],
                                                        pCsIMedium,
                                                        halfLengthCry_CEPA[i],
                                                        vertices_Cry_CEPA[j]);
                Cry_vol_CEPA[j]->SetLineColor(kMagenta);
                Cry_vol_CEPA[j]->SetVisLeaves(kTRUE);
                Cry_vol_CEPA[j]->SetVisibility(kTRUE);
                Cry_vol_CEPA[j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[j] =
                    gGeoManager->MakeArb8(WrapCryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA[j],
                                          pWrappingMedium,
                                          halfLengthCry_CEPA[i] + wrapping_thickness_CEPA,
                                          vertices_Cry_Wrap_CEPA[j]);
                Cry_vol_wrap_CEPA[j]->SetLineColor(kGreen);
                Cry_vol_wrap_CEPA[j]->SetVisLeaves(kTRUE);
                Cry_vol_wrap_CEPA[j]->SetVisibility(kTRUE);
                Cry_vol_wrap_CEPA[j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[j]->AddNode(Cry_vol_CEPA[j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
                Alv_inner_vol_CEPA[i]->AddNode(Cry_vol_wrap_CEPA[j],
                                               0,
                                               new TGeoCombiTrans(cry_position_local_CEPA[j].X(),
                                                                  cry_position_local_CEPA[j].Y(),
                                                                  cry_position_local_CEPA[j].Z(),
                                                                  rotUni));
            }
        }
        if (i == 1) // 6 crystal alveoli
        {
            for (Int_t j = 0; j < 6; j++)
            {
                Cry_vol_CEPA[4 + j] = gGeoManager->MakeArb8(CryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA_6[j],
                                                            pCsIMedium,
                                                            halfLengthCry_CEPA[i],
                                                            vertices_Cry_CEPA[4 + j]);
                Cry_vol_CEPA[4 + j]->SetLineColor(kMagenta);
                Cry_vol_CEPA[4 + j]->SetVisLeaves(kTRUE);
                Cry_vol_CEPA[4 + j]->SetVisibility(kTRUE);
                Cry_vol_CEPA[4 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[4 + j] =
                    gGeoManager->MakeArb8(WrapCryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA_6[j],
                                          pWrappingMedium,
                                          halfLengthCry_CEPA[i] + wrapping_thickness_CEPA,
                                          vertices_Cry_Wrap_CEPA[4 + j]);
                Cry_vol_wrap_CEPA[4 + j]->SetLineColor(kGreen);
                Cry_vol_wrap_CEPA[4 + j]->SetVisLeaves(kTRUE);
                Cry_vol_wrap_CEPA[4 + j]->SetVisibility(kTRUE);
                Cry_vol_wrap_CEPA[4 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[4 + j]->AddNode(Cry_vol_CEPA[4 + j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
                Alv_inner_vol_CEPA[i]->AddNode(Cry_vol_wrap_CEPA[4 + j],
                                               0,
                                               new TGeoCombiTrans(cry_position_local_CEPA[4 + j].X(),
                                                                  cry_position_local_CEPA[4 + j].Y(),
                                                                  cry_position_local_CEPA[4 + j].Z(),
                                                                  rotUni));
            }
        }
        if (i == 2) // 4 crystal alveoli
        {
            for (Int_t j = 0; j < 4; j++)
            {
                Cry_vol_CEPA[10 + j] = gGeoManager->MakeArb8(CryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA[j],
                                                             pCsIMedium,
                                                             halfLengthCry_CEPA[i],
                                                             vertices_Cry_CEPA[10 + j]);
                Cry_vol_CEPA[10 + j]->SetLineColor(kMagenta);
                Cry_vol_CEPA[10 + j]->SetVisLeaves(kTRUE);
                Cry_vol_CEPA[10 + j]->SetVisibility(kTRUE);
                Cry_vol_CEPA[10 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[10 + j] =
                    gGeoManager->MakeArb8(WrapCryGlobalName_CEPA + name_Alv_CEPA[i] + name_Cry_CEPA[j],
                                          pWrappingMedium,
                                          halfLengthCry_CEPA[i] + wrapping_thickness_CEPA,
                                          vertices_Cry_Wrap_CEPA[10 + j]);
                Cry_vol_wrap_CEPA[10 + j]->SetLineColor(kGreen);
                Cry_vol_wrap_CEPA[10 + j]->SetVisLeaves(kTRUE);
                Cry_vol_wrap_CEPA[10 + j]->SetVisibility(kTRUE);
                Cry_vol_wrap_CEPA[10 + j]->SetVisContainers(kTRUE);

                Cry_vol_wrap_CEPA[10 + j]->AddNode(Cry_vol_CEPA[10 + j], 0, new TGeoCombiTrans(0, 0, 0, rotUni));
                Alv_inner_vol_CEPA[i]->AddNode(Cry_vol_wrap_CEPA[10 + j],
                                               0,
                                               new TGeoCombiTrans(cry_position_local_CEPA[10 + j].X(),
                                                                  cry_position_local_CEPA[10 + j].Y(),
                                                                  cry_position_local_CEPA[10 + j].Z(),
                                                                  rotUni));
            }
        }

        // Inner volume center is displaced cf_thickness_CEPA microns along Z
        Alv_vol_CEPA[i]->AddNode(Alv_inner_vol_CEPA[i], 0, new TGeoCombiTrans(0, 0, cf_thickness_CEPA, rotUni));

        rotEle_CEPA[0] = rot_CEPA[2 * i].XX();
        rotEle_CEPA[1] = rot_CEPA[2 * i].XY();
        rotEle_CEPA[2] = rot_CEPA[2 * i].XZ();
        rotEle_CEPA[3] = rot_CEPA[2 * i].YX();
        rotEle_CEPA[4] = rot_CEPA[2 * i].YY();
        rotEle_CEPA[5] = rot_CEPA[2 * i].YZ();
        rotEle_CEPA[6] = rot_CEPA[2 * i].ZX();
        rotEle_CEPA[7] = rot_CEPA[2 * i].ZY();
        rotEle_CEPA[8] = rot_CEPA[2 * i].ZZ();
        rotAlv_CEPA[i]->SetMatrix(rotEle_CEPA);

        for (Int_t j = 0; j < 8; j++)
        { // rotation around Z
            // rotAlvFinal_CEPA[i * 32 + j] = new TGeoRotation((*rotOnZ_CEPA[j]) * (*rotAlv_CEPA[i]));
            // alv_cm_rot_CEPA[2 * i] = (*rotationOnZ_CEPA[j]) * alv_cm_CEPA[2 * i];

            rotAlvFinal_CEPA[i * 8 + j] = new TGeoRotation((*rotOnZ_CEPA[j]) * (*rotAlv_CEPA[i]));
            alv_cm_rot_CEPA[i] = (*rotationOnZ_CEPA[j]) * alv_cm_CEPA[i];

            if (isCrystalInstalled(i + 24, j, installedCrystals)) // alveoli number stars in 24 for CEPA
            {
                pWorld->AddNode(Alv_vol_CEPA[i],
                                j,
                                new TGeoCombiTrans(alv_cm_rot_CEPA[i].X(),
                                                   alv_cm_rot_CEPA[i].Y(),
                                                   alv_cm_rot_CEPA[i].Z(),
                                                   rotAlvFinal_CEPA[i * 8 + j]));
            }
        }
    }
  
    ConstructGeometry(pMedVac,pMed1,pMed2,pMed3,pMed4,pMed5,pMedSi,pMedGold,pMedPolyethylene,pMed41,pMed42,pMed43);

    gGeoMan->CloseGeometry();
    gGeoMan->CheckOverlaps(0.001);
    gGeoMan->PrintOverlaps();
    gGeoMan->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();

    std::cout << "\033[34m Creating geometry:\033[0m "
              << "\033[33m" << geoFileName << " \033[0m" << std::endl;
    std::cout << "Macro finished successfully." << std::endl;
}

Bool_t isCrystalInstalled(Int_t alvType, Int_t alveolusCopy, Int_t instCry[])
{
    // reproduces partially the algorithm of R3BCalifaGeometry::GetCrystalId(const char* volumePath)
    Bool_t found = kFALSE;
    Int_t crystalId = 0;
    Int_t cryType = 1; // first crystal of the alveoli... if not present, alveoli is removed.

    if (alvType == 1)
        crystalId = 1 + alveolusCopy; // first alveoli ring, one crystal per alveolus
    else if (alvType < 20)
        crystalId = 33 + (alvType - 2) * 128 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus
    else if (alvType < 24)
        crystalId = 2337 + (alvType - 20) * 24 + alveolusCopy * 3 + (cryType - 1); // three crystal per alveolus
    else if (alvType == 24)
        crystalId = 2433 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus (2433 to 2464)
    else if (alvType == 25)
        crystalId = 2465 + alveolusCopy * 6 + (cryType - 1); // six crystal per alveolus (2465 to 2512)
    else if (alvType == 26)
        crystalId = 2513 + alveolusCopy * 4 + (cryType - 1); // four crystal per alveolus (2513 to 2544)
    else
    {
        cout << "isCrystalInstalled: Wrong alveolus number ";
        cout << "---- alvType: " << alvType << endl;
        return 0;
    }

    if (crystalId < 1 || crystalId > nbcrystals)
    { // crystalId runs from 1 to nbcrystals
        cout << "isCrystalInstalled: Wrong crystal number ";
        cout << "---- crystalId: " << crystalId << endl;
        return 0;
    }
    for (Int_t i = 0; i < nbcrystals; i++)
    {
        if (crystalId == instCry[i])
            found = kTRUE;
    }
    return found;
}

void ConstructGeometry(TGeoMedium *pMedVac, TGeoMedium *pMed1,
		       TGeoMedium *pMed2, TGeoMedium *pMed3,
		       TGeoMedium *pMed4, TGeoMedium *pMed5,
		       TGeoMedium *pMedSi, TGeoMedium *pMedGold,
		       TGeoMedium *pMedPolyethylene, TGeoMedium *pMed41, TGeoMedium *pMed42, TGeoMedium *pMed43){
  cout << endl;
  cout << "-I- R3BTarget:: ConstructGeometry() "<< endl;
  cout << "-I- R3BTarget Target type: LiH "<< endl;
  cout << endl;
  
  // Shape: World type: TGeoBBox
  
  // Creating enveloppe
  TGeoShape* pVCTube = new TGeoTube("CenterPart", 0., 29.26, 140.2 / 2.0);

  TGeoCombiTrans* t_tube = new TGeoCombiTrans("t_tube", 0., 0., 0., fRefRot);
  t_tube->RegisterYourself();

  TGeoShape* pVCCone = new TGeoCone("Cone", (33.688 + 4.) / 2., 0., 4.6, 0., 7.8);

  TGeoCombiTrans* t_cone = new TGeoCombiTrans("t_cone", 0., 0., 442. / 2.0 + (33.688 + 4.) / 2., fRefRot);
  t_cone->RegisterYourself();

  TGeoCompositeShape* pVCWorld = new TGeoCompositeShape("VCbox", "CenterPart:t_tube + Cone:t_cone");
    
  TGeoVolume* pTargetEnveloppe_log = new TGeoVolume("TargetEnveloppe", pVCWorld, pMedVac);
  pTargetEnveloppe_log->SetVisLeaves(kTRUE);
  TGeoCombiTrans* t0 = new TGeoCombiTrans();
  TGeoCombiTrans* pGlobal = GetGlobalPosition(t0);
  //pWorld->AddNode(pTargetEnveloppe_log, 0, pGlobal);

  Double_t dx, dy, dz;
  Double_t thx, thy, thz;
  Double_t phx, phy, phz;

  Double_t offsetZ = -0.50; // offset with respect to LH2 target cell, entrance at (0,0,0)

  // Mylar windows for LH2 target
  Double_t thinMylar1 = 0.012 / 2.0;            // 120micra entrance window
  Double_t thinMylar2 = 0.018 / 2.0;            // 180micra exit window
  Double_t LH2targetR = 1.5 - thinMylar2 * 2.0; // Radius of 1.5 cm - mylar cylindre
  //Double_t LH2targetL = 1.5;                  // Length 1.5cm
  Double_t LH2targetL = 5.;                     // Length 5.cm

  //-------- Chamber -------------------------------------------------------------------------
  Double_t parChamberBarrel[3] = { 26.1, 26.4, (46.3 / 2.) };
  TGeoTube* pChamberBarrel = new TGeoTube(parChamberBarrel);
  TGeoVolume* pChamberBarrelLog = new TGeoVolume("ChamberBarrelLog", pChamberBarrel, pMed2);
  pChamberBarrelLog->SetVisLeaves(kTRUE);
  pChamberBarrelLog->SetLineColor(34);

  dx = 0.000000;
  dy = 0.000000;
  dz = -42.1250 + 1.7 + 46.3 / 2. + offsetZ;
  TGeoRotation* rot0 = new TGeoRotation();
  TGeoCombiTrans* pMatrix169 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberBarrelLog, 0, pMatrix169);

  // chamber Inner
  Double_t parChamberBarrelIn[3] = { 9.625, (9.625 + 2.2), (37.57 / 4.) };
  TGeoTube* pChamberBarrelIn = new TGeoTube(parChamberBarrelIn);
  TGeoVolume* pChamberBarrelInLog = new TGeoVolume("ChamberBarrelInLog", pChamberBarrelIn, pMed2);
  pChamberBarrelInLog->SetVisLeaves(kTRUE);
  pChamberBarrelInLog->SetLineColor(34);

  dx = 0.000000;
  dy = 0.000000;
  dz = -37.54000 - 37.57 / 4. + offsetZ;
  TGeoCombiTrans* pMatrix170 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberBarrelInLog, 0, pMatrix170);

  // Inner support
  Double_t parChamberIn[3] = { 0., 0.6, ((37.57 - 0.9) / 2.) };
  TGeoTube* pChamberIn = new TGeoTube(parChamberIn);
  TGeoVolume* pChamberInLog = new TGeoVolume("ChamberInLog", pChamberIn, pMed2);
  pChamberInLog->SetVisLeaves(kTRUE);
  pChamberInLog->SetLineColor(34);

  dx = 9.6 - 0.6;
  dy = 0.000000;
  dz = -37.54000 + 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix100 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberInLog, 0, pMatrix100);

  dx = -(9.6 - 0.6);
  dy = 0.000000;
  dz = -37.54000 + 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix101 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberInLog, 1, pMatrix101);

  dy = -(9.6 - 0.6);
  dx = 0.000000;
  dz = -37.54000 + 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix102 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberInLog, 2, pMatrix102);

  dy = (9.6 - 0.6);
  dx = 0.000000;
  dz = -37.54000 + 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix103 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberInLog, 3, pMatrix103);

  Double_t parChamberIn2[3] = { 7.6, 9.625, (0.9 / 2.) };
  TGeoTube* pChamberIn2 = new TGeoTube(parChamberIn2);
  TGeoVolume* pChamberIn2Log = new TGeoVolume("ChamberIn2Log", pChamberIn2, pMed2);
  pChamberIn2Log->SetVisLeaves(kTRUE);
  pChamberIn2Log->SetLineColor(17);

  dx = 0.000000;
  dy = 0.000000;
  dz = -17.8550 - 0.9 - 37.57 + 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix104 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberIn2Log, 3, pMatrix104);

  Double_t parChamberEndPlate[3] = { (9.625 + 2.2), 29.25, 1.7 };
  TGeoTube* pChamberEndPlate = new TGeoTube(parChamberEndPlate);
  TGeoVolume* pChamberEndPlateLog = new TGeoVolume("ChamberEndPlateLog", pChamberEndPlate, pMed3);
  pChamberEndPlateLog->SetVisLeaves(kTRUE);
  pChamberEndPlateLog->SetLineColor(34);

  dx = 0.000000;
  dy = 0.000000;
  dz = -42.1250 + offsetZ;
  TGeoCombiTrans* pMatrix174 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberEndPlateLog, 0, pMatrix174);

  Double_t parChamberEndPlateIn[3] = { (3.8), (9.625 + 2.2), 0.45 };
  TGeoTube* pChamberEndPlateIn = new TGeoTube(parChamberEndPlateIn);
  TGeoVolume* pChamberEndPlateInLog = new TGeoVolume("ChamberEndPlateInLog", pChamberEndPlateIn, pMed2);
  pChamberEndPlateInLog->SetVisLeaves(kTRUE);
  pChamberEndPlateInLog->SetLineColor(34);

  dx = 0.000000;
  dy = 0.000000;
  dz = -17.8550 - 0.45 + offsetZ;
  TGeoCombiTrans* pMatrix175 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pChamberEndPlateInLog, 0, pMatrix175);

  Double_t parTargetEndPlateIn[3] = { (LH2targetR + thinMylar2 * 2.0), 2.44, 0.75 };
  TGeoTube* pTargetEndPlateIn = new TGeoTube(parTargetEndPlateIn);
  TGeoVolume* pTargetEndPlateInLog = new TGeoVolume("TargetEndPlateLog", pTargetEndPlateIn, pMed1);
  pTargetEndPlateInLog->SetVisLeaves(kTRUE);
  pTargetEndPlateInLog->SetLineColor(46);

  dx = 0.000000;
  dy = 0.000000;
  dz = 0.75000 + offsetZ;
  TGeoCombiTrans* pMatrix176 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pTargetEndPlateInLog, 0, pMatrix176);

  TGeoSphere* SphereChamber = new TGeoSphere("SphereChamber", 26.1, 26.4, 15.26, 90., 0., 360.);

  dx = 0.000000;
  dy = 0.000000;
  dz = 6.775000 - 0.9;
  TGeoCombiTrans* pMatrix177 = new TGeoCombiTrans("", dx, dy, dz, rot0);

  Double_t parTargetEndPlateOut[3] = { 7.0, 10., 1. };
  TGeoTube* pTargetEndPlateOut = new TGeoTube(parTargetEndPlateOut);
  TGeoVolume* pTargetEndPlateOutLog = new TGeoVolume("TargetEndPlateOutLog", pTargetEndPlateOut, pMed2);
  pTargetEndPlateOutLog->SetVisLeaves(kTRUE);
  pTargetEndPlateOutLog->SetLineColor(34);

  dx = 0.000000;
  dy = 0.000000;
  dz = 6.775000 - 0.9 + 28.812028 + 33.688 + 1. + offsetZ;
  TGeoCombiTrans* pMatrix179 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pTargetEndPlateOutLog, 0, pMatrix179);

  Double_t parRingOut[5] = { (4.301 + 2.348 + 0.1), 2.048 + 0.3, 2.348 + 0.3, 0., 360. };
  TGeoTorus* pRingOut = new TGeoTorus(parRingOut);

  Double_t parCone[5] = { (1.9), 4.601, 10., 4.601, 10. };
  TGeoCone* pCone = new TGeoCone(parCone);

  dx = 0.000000;
  dy = 0.000000;
  dz = 1.;
  TGeoCombiTrans* pMatrixSub = new TGeoCombiTrans("", dx, dy, dz, rot0);

  TGeoSubtraction* pSub1 = new TGeoSubtraction(pRingOut, pCone, 0, pMatrixSub);

  TGeoSphere* Sphere = new TGeoSphere("Sphere", 26.4, 30., 15.2, 25., 0., 360.);

  TGeoShape* pROutb = new TGeoCompositeShape("pROutb", pSub1);

  dx = 0.000000;
  dy = 0.000000;
  dz = -26.4 - 1.5;
  TGeoCombiTrans* pMatrixSub1 = new TGeoCombiTrans("", dx, dy, dz, rot0);

  TGeoSubtraction* pSub2 = new TGeoSubtraction(pROutb, Sphere, 0, pMatrixSub1);

  TGeoShape* pRingOutb1 = new TGeoCompositeShape("pRingOutb1", pSub2);

  Double_t parChamberCone[5] = { (33.688 / 2.), 4.301, 4.601, 7.0, 7.3 };
  TGeoCone* pChamberCone = new TGeoCone(parChamberCone);

  dx = 0.000000;
  dy = 0.000000;
  dz = 6.775000 - 0.9 + 28.812028 + 33.688 / 2.0;
  TGeoCombiTrans* pMatrix178 = new TGeoCombiTrans("", dx, dy, dz, rot0);

  dx = 0.0000;
  dy = 0.0000;
  dz = 6.775000 - 0.9 + 28.812028 - 2.348 / 2. + 0.15;
  TGeoCombiTrans* pMatrix180 = new TGeoCombiTrans("", dx, dy, dz, rot0);

  TGeoUnion* pBoolNode0 = new TGeoUnion(pChamberCone, pRingOutb1, pMatrix178, pMatrix180);

  TGeoShape* pRingOutb = new TGeoCompositeShape("pNode0", pBoolNode0);

  TGeoUnion* pBoolNode1 = new TGeoUnion(pRingOutb, SphereChamber, 0, pMatrix177);

  TGeoShape* pFrontChamber = new TGeoCompositeShape("pFrontChamber", pBoolNode1);

  TGeoVolume* pFrontChamber_log = new TGeoVolume("FrontChamber", pFrontChamber, pMed2);
  pFrontChamber_log->SetVisLeaves(kTRUE);
  pFrontChamber_log->SetLineColor(34);

  dx = 0.0000;
  dy = 0.0000;
  dz = offsetZ;
  TGeoCombiTrans* pMatrix220 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pFrontChamber_log, 0, pMatrix220);

  Double_t parCoolRingIn[5] = { 2., 0.2, 0.35, 0., 90. };
  TGeoTorus* pCoolRingIn = new TGeoTorus(parCoolRingIn);
  TGeoVolume* pCoolRingIn_log = new TGeoVolume("CoolRingIn", pCoolRingIn, pMed1);
  pCoolRingIn_log->SetVisLeaves(kTRUE);
  pCoolRingIn_log->SetLineColor(46);

  dx = 0.000000;
  dy = 2.440000;
  dz = 0.75 - 2. + offsetZ;
  TGeoRotation* rot1 = new TGeoRotation();
  rot1->RotateY(-90.0);
  TGeoCombiTrans* pMatrix5 = new TGeoCombiTrans("", dx, dy, dz, rot1);
  pWorld->AddNode(pCoolRingIn_log, 0, pMatrix5);

  Double_t parCoolRingIn2[5] = { 2., 0.2, 0.35, 270., 90. };
  TGeoTorus* pCoolRingIn2 = new TGeoTorus(parCoolRingIn2);
  TGeoVolume* pCoolRingIn2_log = new TGeoVolume("CoolRingIn2", pCoolRingIn2, pMed1);
  pCoolRingIn2_log->SetVisLeaves(kTRUE);
  pCoolRingIn2_log->SetLineColor(46);

  dx = 0.000000;
  dy = -2.44000;
  dz = 0.75 - 2. + offsetZ;
  TGeoCombiTrans* pMatrix6 = new TGeoCombiTrans("", dx, dy, dz, rot1);
  pWorld->AddNode(pCoolRingIn2_log, 0, pMatrix6);

  Double_t parCoolPipeIn[3] = { 0.2, 0.35, ((17.855 - 1.25) / 2.) };
  TGeoTube* pCoolPipeIn = new TGeoTube(parCoolPipeIn);
  TGeoVolume* pCoolPipeIn_log = new TGeoVolume("CoolPipeIn", pCoolPipeIn, pMed1);
  pCoolPipeIn_log->SetVisLeaves(kTRUE);
  pCoolPipeIn_log->SetLineColor(46);

  dx = 0.000000;
  dy = -2.44000 - 2.;
  dz = -(17.855 - 1.25) / 2. - 1.25 + offsetZ;
  TGeoCombiTrans* pMatrix7 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pCoolPipeIn_log, 0, pMatrix7);

  dy = 2.44000 + 2.;
  TGeoCombiTrans* pMatrix8 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pCoolPipeIn_log, 1, pMatrix8);

  // LH2 target --------------------------------------------------------------
  Double_t parTargetH2[3] = { 0., LH2targetR, (LH2targetL / 2.0) };
  TGeoTube* pTargetH2 = new TGeoTube(parTargetH2);
  TGeoVolume* pTargetH2Log = new TGeoVolume("TargetH2Log", pTargetH2, pMed5);
  pTargetH2Log->SetVisLeaves(kTRUE);
  pTargetH2Log->SetLineColor(1);

  dx = 0.000000;
  dy = 0.000000;
  dz = LH2targetL / 2.0;
  TGeoCombiTrans* pMatrix9 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pTargetH2Log, 0, pMatrix9);

  // Mylar windows for LH2 target
  Double_t parMylar1[3] = { 0., LH2targetR + thinMylar2 * 2., thinMylar1 };
  TGeoTube* pMylar1 = new TGeoTube(parMylar1);
  TGeoVolume* pMylar1Log = new TGeoVolume("MylarWindowEntranceLog", pMylar1, pMed4);
  pMylar1Log->SetVisLeaves(kTRUE);
  pMylar1Log->SetLineColor(12);

  dx = 0.000000;
  dy = 0.000000;
  dz = -thinMylar1;
  TGeoCombiTrans* pMatrix10 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pMylar1Log, 0, pMatrix10);

  Double_t parMylar2[3] = { 0., LH2targetR + thinMylar2 * 2., thinMylar2 };
  TGeoTube* pMylar2 = new TGeoTube(parMylar2);
  TGeoVolume* pMylar2Log = new TGeoVolume("MylarWindowExitLog", pMylar2, pMed4);
  pMylar2Log->SetVisLeaves(kTRUE);
  pMylar2Log->SetLineColor(12);
  dx = 0.000000;
  dy = 0.000000;
  dz = LH2targetL + thinMylar2;
  TGeoCombiTrans* pMatrix11 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pMylar2Log, 0, pMatrix11);

  Double_t parMylar3[3] = { LH2targetR, LH2targetR + thinMylar2 * 2., (LH2targetL / 2.0) };
  TGeoTube* pMylar3 = new TGeoTube(parMylar3);
  TGeoVolume* pMylar3Log = new TGeoVolume("MylarCylindreLog", pMylar3, pMed4);
  pMylar3Log->SetVisLeaves(kTRUE);
  pMylar3Log->SetLineColor(12);
  dx = 0.000000;
  dy = 0.000000;
  dz = LH2targetL / 2.0;
  TGeoCombiTrans* pMatrix12 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pMylar3Log, 0, pMatrix12);

  //Insulation of the LH2 target
  Double_t InsulationThickness = 0.0050; //10*5um = 50um
  Double_t parMylar4[3] = { 0., LH2targetR + thinMylar2 * 2. + InsulationThickness, InsulationThickness/2.};
  TGeoTube* pMylar4 = new TGeoTube(parMylar4);
  TGeoVolume* pMylar4Log = new TGeoVolume("InsulationMylarWindowExitLog", pMylar4, pMed4);
  pMylar4Log->SetVisLeaves(kTRUE);
  pMylar4Log->SetLineColor(12);
  dx = 0.000000;
  dy = 0.000000;
  dz = LH2targetL + thinMylar2 * 2. + InsulationThickness/2.;
  TGeoCombiTrans* pMatrix14 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pMylar4Log, 0, pMatrix14);

  Double_t parMylar5[3] = { LH2targetR + thinMylar2 * 2., LH2targetR + thinMylar2 * 2. + InsulationThickness, thinMylar2 + (LH2targetL / 2.0) - 1./2.};
  TGeoTube* pMylar5 = new TGeoTube(parMylar5);
  TGeoVolume* pMylar5Log = new TGeoVolume("InsulationMylarCylindreLog", pMylar5, pMed4);
  pMylar5Log->SetVisLeaves(kTRUE);
  pMylar5Log->SetLineColor(12);
  dx = 0.000000;
  dy = 0.000000;
  //dz = (LH2targetL) / 2.0 + thinMylar2;
  dz = thinMylar2 + (LH2targetL / 2.0) - 1./2. + 1.;
  TGeoCombiTrans* pMatrix15 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pMylar5Log, 0, pMatrix15);
  
  // Exit mylar window for the vacuum chamber
  Double_t parWin[3] = { 0., 7., 0.0100 / 2.0 }; // 100 micras
  TGeoTube* pWin1 = new TGeoTube(parWin);
  TGeoVolume* pWin1Log = new TGeoVolume("ExitWinLog", pWin1, pMed4);
  pWin1Log->SetVisLeaves(kTRUE);
 // pWin1Log->SetLineColor(2);

  dx = 0.000000;
  dy = 0.000000;
  dz = 6.775000 - 0.9 + 28.812028 + 33.688 + 1. + offsetZ;
  TGeoCombiTrans* pMatrix13 = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pWin1Log, 0, pMatrix13);
  
  //FOOT Definition
  //Al support
  // TGeoShape* support_foot_box = new TGeoBBox("support_foot_box",
  // 					 12., 
  // 					 0.4, 
  // 					 20.);
  
  // TGeoVolume *supportfootvol = new TGeoVolume("supportfootvol", support_foot_box , pMed2);
  // supportfootvol->SetVisLeaves(kTRUE);
  // supportfootvol->SetLineColor(kBlue);

  // dx = 0.000000;
  // dy = -15.0000;
  // dz = 0.000000;
  // TGeoCombiTrans* pMatrixSF = new TGeoCombiTrans("", dx, dy, dz, rot0);
  // pWorld->AddNode(supportfootvol, 0, pMatrixSF);
  
  //Test
  /*
  Double_t parTest[3] = { 0.,2.,4};
  TGeoTube* pTest = new TGeoTube(parTest);
  TGeoVolume* pTestLog = new TGeoVolume("TestLog", pTest, pMed4);
  pTestLog->SetVisLeaves(kTRUE);
  pTestLog->SetLineColor(12);

  dx = 0.000000;
  dy = 0.000000;
  dz = 15.000000;
  TGeoCombiTrans* pMatrixTest = new TGeoCombiTrans("", dx, dy, dz, rot0);
  pWorld->AddNode(pTestLog, 0, pMatrixTest);
  */
  
  //Positions
  // Double_t PosX1 = 5.2; //X position first layer
  // Double_t PosZ1 = 2.5+3.5; //Z position first layer
  // Double_t Theta1 = 65.; //65 //56.6 //angle first layer
  // Double_t Gap1 = 0.5;
  
  // Double_t PosX2 = 7.; //X position second layer
  // Double_t PosZ2 = 11.; //Z position second layer
  // Double_t Theta2 = 49; //49 //angle second layer
  // Double_t Gap2 = 0.5;

  //Parallel double FOOTS
  // Double_t PosX1 = 5.2; //X position first layer
  // Double_t PosZ1 = 2.5; //Z position first layer
  // Double_t Theta1 = TMath::RadToDeg()*TMath::ATan(PosX1/PosZ1); //65 //56.6 //angle first layer
  // Double_t Gap1 = 0.5;
  
  // Double_t PosX2 = 7.; //X position second layer
  // Double_t PosZ2 = 11.; //Z position second layer
  // Double_t Theta2 = TMath::RadToDeg()*TMath::ATan(PosX2/PosZ2); //49 //angle second layer
  // Double_t Gap2 = 0.5;

  //FOOT Aligned on a certain Angle Trajectory
  Double_t alpha = 18.; //cover 20° trajectories
  Double_t ATheta = 65.;
  Double_t Theta1 = 65.;
  Double_t Theta2 = 65.;
  double shift = 0.;
  double shiftX = 0.5;
  
  //Double_t PosZ1 = 6.; //X position first layer
  //Double_t Theta1 = 60; //56.6 //angle first layer
  Double_t Theta1_comp = 90. - Theta1;
  double fdz1 = 8.3; //5.2cm
  Double_t GapAF = 1.2;
  Double_t PosZ1 = fdz1 - (5.0 -4.5) * TMath::Cos(TMath::DegToRad()*Theta1) + GapAF/TMath::Cos(TMath::DegToRad()*(Theta1-alpha))*TMath::Cos(TMath::DegToRad()*alpha);
  
  Double_t PosX1 = 5.*TMath::Sin(TMath::DegToRad()*Theta1_comp)+TMath::Tan(TMath::DegToRad()*alpha)*((PosZ1-2.5)+5*TMath::Cos(TMath::DegToRad()*Theta1_comp)); //Z position first layer
  Double_t Gap1 = 0.82; //0.5+2x0.16(PCBs)
  Double_t PosZ1b = (1./(1.+TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta1_comp)))*(PosZ1+PosX1/TMath::Tan(TMath::DegToRad()*Theta1_comp)- 5.*TMath::Sin(TMath::DegToRad()*Theta1_comp)/TMath::Tan(TMath::DegToRad()*Theta1_comp)+Gap1/TMath::Sin(TMath::DegToRad()*Theta1_comp)- 5.*TMath::Cos(TMath::DegToRad()*Theta1_comp)*TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta1_comp)+2.5*TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta1_comp));
  
 // PosZ1b = PosZ1 + Gap1/TMath::Cos(TMath::DegToRad()*Theta1)/(1.+TMath::Tan(TMath::DegToRad()*alpha)*TMath::Tan(TMath::DegToRad()*Theta1));
  PosZ1b = PosZ1 + Gap1/TMath::Cos(TMath::DegToRad()*(Theta1-alpha))*TMath::Cos(TMath::DegToRad()*alpha);
 
  //TMath::ATan((PosX1b-5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+5*TMath::Cos(TMath::DegToRad()*Theta1_comp)+3.5(double) 
  Double_t PosX1b = 5.*TMath::Sin(TMath::DegToRad()*Theta1_comp)+TMath::Tan(TMath::DegToRad()*alpha)*((PosZ1b-2.5)+5*TMath::Cos(TMath::DegToRad()*Theta1_comp));

  Double_t GapF = 2.0;  
  Double_t PosZ2 = GapF/TMath::Cos(TMath::DegToRad()*(Theta1-alpha))*TMath::Cos(TMath::DegToRad()*alpha) + PosZ1;
 // Double_t PosZ2 = 10.5; //X position second layer
  //Double_t Theta2 = 60; //56.6 //angle second layer
  Double_t Theta2_comp = 90. - Theta2;
  Double_t PosX2 = 5.*TMath::Sin(TMath::DegToRad()*Theta2_comp)+TMath::Tan(TMath::DegToRad()*alpha)*((PosZ2-2.5)+5*TMath::Cos(TMath::DegToRad()*Theta2_comp)); //Z position second layer
  Double_t Gap2 = 0.82; //0.5+2x0.16(PCBs)
  Double_t PosZ2b = (1./(1.+TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta2_comp)))*(PosZ2+PosX2/TMath::Tan(TMath::DegToRad()*Theta2_comp)-5.*TMath::Sin(TMath::DegToRad()*Theta2_comp)/TMath::Tan(TMath::DegToRad()*Theta2_comp)+Gap2/TMath::Sin(TMath::DegToRad()*Theta2_comp)-5.*TMath::Cos(TMath::DegToRad()*Theta2_comp)*TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta2_comp)+2.5*TMath::Tan(TMath::DegToRad()*alpha)/TMath::Tan(TMath::DegToRad()*Theta2_comp));
  Double_t PosX2b = 5.*TMath::Sin(TMath::DegToRad()*Theta2_comp)+TMath::Tan(TMath::DegToRad()*alpha)*((PosZ2b-2.5)+5*TMath::Cos(TMath::DegToRad()*Theta2_comp));
  
  
  Double_t PosZ3 = 30.; //Z position beam tracker
  Double_t Gap3 = 0.5;
  
  cout << "/*****************************************/" << endl;
  cout << "/*************FOOT Layer 1****************/" << endl;
  cout << "Pos X1a : " << PosX1 << " cm" << endl;
  cout << "Pos Z1a : " << PosZ1 << " cm" << endl;
  cout << "Dist 1a : " << TMath::Sqrt(PosX1*PosX1+(PosZ1+shift-2.5)*(PosZ1+shift-2.5)) << " cm" << endl;
  cout << "Pos X1b : " << PosX1b << " cm" << endl;
  cout << "Pos Z1b : " << PosZ1b << " cm" << endl;
  cout << "Dist 1b : " << TMath::Sqrt(PosX1b*PosX1b+(PosZ1b+shift-2.5)*(PosZ1b+shift-2.5)) << " cm" << endl;
  cout << "Theta 1 : " << Theta1 << " deg" << endl;
  cout << "/*************FOOT Layer 2****************/" << endl;
  cout << "Pos X2a : " << PosX2 << " cm" << endl;
  cout << "Pos Z2a : " << PosZ2 << " cm" << endl;
  cout << "Dist 2a : " << TMath::Sqrt(PosX2*PosX2+(PosZ2+shift-2.5)*(PosZ2+shift-2.5)) << " cm" << endl;
  cout << "Pos X2b : " << PosX2b << " cm" << endl;
  cout << "Pos Z2b : " << PosZ2b << " cm" << endl;
  cout << "Dist 2b : " << TMath::Sqrt(PosX2b*PosX2b+(PosZ2b+shift-2.5)*(PosZ2b+shift-2.5)) << " cm" << endl;
  cout << "Theta 2 : " << Theta2 << " deg" << endl;
  cout << "Theta exit : " << TMath::ATan((PosX2b - 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() <<" deg" << endl;
  cout << "Theta Inner : " << TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << " deg" << endl;
  cout << "/*****************************************/" << endl;  
  
  Double_t tx,ty,tz;
  
  TGeoRotation * zeroRot = new TGeoRotation; //zero rotation
  TGeoCombiTrans * tZero = new TGeoCombiTrans("tZero", 0., 0., 0., zeroRot);
  tZero->RegisterYourself();
  
  //------------------ PCBs------------------------------------------
  TGeoShape* padle_h_box1 = new TGeoBBox("padle_h_box1",
					 7.25, 
					 9.1, 
					 0.0750);
  TGeoShape* padle_h_box2 = new TGeoBBox("padle_h_box2",
					 5., 
					 5., 
					 0.0750);
  TGeoShape* padle_h_box3 = new TGeoBBox("padle_h_box3",
					 9.1, 
					 7.25, 
					 0.0750);
  
  //Shift the Si position
  TGeoCombiTrans * tShift = new TGeoCombiTrans("tShift", 0.0, -4.1 + 2.62, 0.075, zeroRot);
  tShift->RegisterYourself();

  TGeoCombiTrans * tShift2 = new TGeoCombiTrans("tShift2", 0.0, 1.48, 0.075, zeroRot);
  tShift2->RegisterYourself();

  TGeoCombiTrans * tShift3 = new TGeoCombiTrans("tShift3", -1.48, 0.0, 0.075, zeroRot);
  tShift3->RegisterYourself();

  TGeoCombiTrans * tShift4 = new TGeoCombiTrans("tShift4", 1.48, 0.0, 0.075, zeroRot);
  tShift4->RegisterYourself();
  
  // Create a composite shape
  TGeoCompositeShape *PCB = new TGeoCompositeShape("diffbox", "padle_h_box1:tShift - padle_h_box2");
  TGeoCompositeShape *PCB2 = new TGeoCompositeShape("diffbox2", "padle_h_box1:tShift2 - padle_h_box2");
  TGeoCompositeShape *PCB3 = new TGeoCompositeShape("diffbox3", "padle_h_box3:tShift3 - padle_h_box2");
  TGeoCompositeShape *PCB4 = new TGeoCompositeShape("diffbox4", "padle_h_box3:tShift4 - padle_h_box2");
  TGeoVolume *PCBvol1 = new TGeoVolume("PCB", PCB, pMedPolyethylene);
  PCBvol1->SetVisLeaves(kTRUE);
  PCBvol1->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol2 = new TGeoVolume("PCB2", PCB2, pMedPolyethylene);
  PCBvol2->SetVisLeaves(kTRUE);
  PCBvol2->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol3 = new TGeoVolume("PCB3", PCB3, pMedPolyethylene);
  PCBvol3->SetVisLeaves(kTRUE);
  PCBvol3->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol4 = new TGeoVolume("PCB4", PCB4, pMedPolyethylene);
  PCBvol4->SetVisLeaves(kTRUE);
  PCBvol4->SetLineColor(kGreen+2);
  
  //**************************************************************//
  //*********************   Si Sensors ******************************//
  //***************************************************************/
  
  // Si Shape & volume: TraBox type: TGeoBBox
  dx = 5.; //0.5*14.50000; //5.
  dy = 5.; //0.5*(10.00000+2.*2.62); //5.
  dz = 0.0075000;
  // Volume: TraLog
  TGeoVolume *TraLog = gGeoManager->MakeBox("TraLog",pMedSi,dx,dy,dz);
  TraLog->SetVisLeaves(kTRUE);
  TraLog->SetLineColor(33);
  
  //TRANSFORMATION MATRICES
  dx = - PosX1;
  dy = 0.000000;
  dz = PosZ1+shift;
  // Rotation:
  thx = 33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *      pMatrix203 = new TGeoRotation();
  pMatrix203->SetAngles(90,-Theta1,0);
  TGeoCombiTrans*  pMatrix202 = new TGeoCombiTrans("", dx,dy,dz,pMatrix203);
  
  
 // TGeoRotation *      alpMatrix203 = new TGeoRotation();
//  alpMatrix203->SetAngles(90,-165,0);
//  TGeoCombiTrans*  alpMatrix202 = new TGeoCombiTrans("", dx+1,dy,dz,alpMatrix203);

  //TRANSFORMATION MATRICES
  // dx = - PosX1*(TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1)+Gap1)/TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1); 
  // dy = 0.000000;
  // dz = (TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1)+Gap1)*PosZ1/TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1); 
  dx = - PosX1b;
  dy = 0.;
  dz = PosZ1b+shift;
  // Rotation:
  thx = 33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *      pMatrix229 = new TGeoRotation();
  pMatrix229->SetAngles(90,-Theta1,0);
  TGeoCombiTrans*  pMatrix228 = new TGeoCombiTrans("", dx,dy,dz,pMatrix229); 
  
  //Combi transformation:
  dx = -PosX2;
  dy = 0.000000;
  dz = PosZ2 + shift;
  // Rotation:
  thx = 41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix217 = new TGeoRotation();
  pMatrix217->SetAngles(90,-Theta2,0);
  TGeoCombiTrans*   pMatrix216 = new TGeoCombiTrans("", dx,dy,dz,pMatrix217);

  //Combi transformation:
  // dx = - PosX2*(TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2)+Gap2)/TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2); 
  // dy = 0.000000;
  // dz = (TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2)+Gap2)*PosZ2/TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2); 
  dx = -PosX2b;
  dy = 0.000000;
  dz = PosZ2b + shift;
  // Rotation:
  thx = 41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix225 = new TGeoRotation();
  pMatrix225->SetAngles(90,-Theta2,0);
  TGeoCombiTrans*   pMatrix224 = new TGeoCombiTrans("", dx,dy,dz,pMatrix225); 
	
  //Combi transformation:
  dx = PosX1;
  dy = 0.000000;
  dz = PosZ1 +shift;
  // Rotation:
  thx = -33.400000;   phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix205 = new TGeoRotation();
  pMatrix205->SetAngles(90,Theta1,0); //56.6
  TGeoCombiTrans*   pMatrix204 = new TGeoCombiTrans("", dx,dy,dz,pMatrix205);


  TGeoRotation *       alpMatrix205 = new TGeoRotation();
  alpMatrix205->SetAngles(90,165,0);  
  TGeoCombiTrans*   alpMatrix204 = new TGeoCombiTrans("", dx-1,dy,dz,alpMatrix205);

  //Combi transformation:
  // dx = PosX1*(TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1)+Gap1)/TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1); 
  // dy = 0.000000;
  // dz = (TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1)+Gap1)*PosZ1/TMath::Sqrt(PosX1*PosX1+PosZ1*PosZ1); 
  dx = PosX1b;
  dy = 0.;
  dz = PosZ1b + shift;
  // Rotation:
  thx = -33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix227 = new TGeoRotation();
  pMatrix227->SetAngles(90,Theta1,0);
  TGeoCombiTrans*   pMatrix226 = new TGeoCombiTrans("", dx,dy,dz,pMatrix227); 

  //Combi transformation:
  dx = PosX2;
  dy = 0.000000;
  dz = PosZ2 + shift;
  // Rotation:
  thx = -41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix219 = new TGeoRotation();
  pMatrix219->SetAngles(90,Theta2,0);
  TGeoCombiTrans*   pMatrix218 = new TGeoCombiTrans("", dx,dy,dz,pMatrix219);

  //Combi transformation:
  // dx = PosX2*(TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2)+Gap2)/TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2); 
  // dy = 0.000000;
  // dz = (TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2)+Gap2)*PosZ2/TMath::Sqrt(PosX2*PosX2+PosZ2*PosZ2); 
  dx = PosX2b;
  dy = 0.000000;
  dz = PosZ2b + shift;
  // Rotation:
  thx = -41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix223 = new TGeoRotation();
  pMatrix223->SetAngles(90,Theta2,0);
  TGeoCombiTrans*   pMatrix222 = new TGeoCombiTrans("", dx,dy,dz,pMatrix223); 
	
  //Combi transformation:
  dx = 0.00000;
  dy = 0.00000;
  dz = PosZ3 + 0.;
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix215 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix214 = new TGeoCombiTrans("", dx,dy,dz,pMatrix215); //Beam Tracker 1

  //Combi transformation:
  dx = 0.00000;
  dy = 0.00000;
  dz = PosZ3+Gap3 + 0.;
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix221 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix2201 = new TGeoCombiTrans("", dx,dy,dz,pMatrix221); //Beam Tracker 2
	
  /************ Assembling everything together ****************/
  TGeoVolume *aTra = new TGeoVolumeAssembly("ATRA");
    
  //PCB
  aTra->AddNode(PCBvol1,11, pMatrix204);
  aTra->AddNode(PCBvol2,12, pMatrix202);
  aTra->AddNode(PCBvol2,13, pMatrix216);
  aTra->AddNode(PCBvol1,14, pMatrix218);
  //aTra->AddNode(PCBvol1,15, pMatrix214); //Tracker
  aTra->AddNode(PCBvol3,16, pMatrix226);
  aTra->AddNode(PCBvol3,17, pMatrix228);
  aTra->AddNode(PCBvol3,18, pMatrix224);
  aTra->AddNode(PCBvol3,19, pMatrix222);
  //aTra->AddNode(PCBvol3,20, pMatrix2201); //Tracker
  
  //Si sensors
  aTra->AddNode(TraLog,1, pMatrix204);
  aTra->AddNode(TraLog,2, pMatrix202);
  aTra->AddNode(TraLog,3, pMatrix216);
  aTra->AddNode(TraLog,4, pMatrix218);
  //aTra->AddNode(TraLog,5, pMatrix214); //Tracker
  aTra->AddNode(TraLog,6, pMatrix226);
  aTra->AddNode(TraLog,7, pMatrix228);
  aTra->AddNode(TraLog,8, pMatrix224);
  aTra->AddNode(TraLog,9, pMatrix222);
  //aTra->AddNode(TraLog,10, pMatrix2201); //Tracker
  
  TGeoRotation *rotg = new TGeoRotation();
  rotg->RotateX(0.);
  rotg->RotateY(0.);
  rotg->RotateZ(0.);
  dx=tx=0.0;
  dy=ty=0.0;
  dz=tz=0.0;
  
  TGeoCombiTrans *t_zero = new TGeoCombiTrans("t_zero");

  TGeoCombiTrans *t00 = new TGeoCombiTrans(tx,ty,tz,rotg);

  pWorld->AddNode(aTra, 1, GetGlobalPosition(t00));


  
  TGeoVolume *alpd = new TGeoVolumeAssembly("ALPD");
  
    // ALPIDE two arm configuration

    float thAlpide = 0.0050; // 50 um

    // Frame
    float Alframelength = 13.4; // cm
    float AlframeWd = 0.4;      // cm
    float AlRingWd = 0.2;       // cm

    auto frame = new TGeoBBox("Frame", Alframelength / 2., Alframelength / 2., AlframeWd / 2.);
    auto frameinner = new TGeoBBox("Frameinner", 10.6 / 2., 10.6 / 2., AlframeWd);
    TGeoSubtraction* rframe = new TGeoSubtraction(frame, frameinner, 0, 0);
    TGeoShape* rframe2 = new TGeoCompositeShape("rframe", rframe);
    TGeoVolume* realframe = new TGeoVolume("AlFrame", rframe2, pMed2);
    realframe->SetVisLeaves(kTRUE);
    realframe->SetLineColor(kBlue);

    auto ring = new TGeoBBox("Ring", Alframelength / 2., Alframelength / 2., AlRingWd / 2.);
    TGeoSubtraction* rring = new TGeoSubtraction(ring, frameinner, 0, 0);
    TGeoShape* rring2 = new TGeoCompositeShape("rring", rring);
    TGeoVolume* realring = new TGeoVolume("AlRing", rring2, pMed2);
    realring->SetVisLeaves(kTRUE);
    realring->SetLineColor(kGreen);
    TGeoVolume* realring1 = new TGeoVolume("AlRing", rring2, pMed2);
   
    
    

 /*-------------------------------------------------------ALPIDE LEFT ARM--------------------------------------------------*/ 
     TGeoVolume* Plane0 = new TGeoVolumeAssembly("APlane0");
  
  TGeoVolume *NTraLog[36]; 
  Char_t buffer[126];
  for(Int_t n = 0; n < 36 ; n++){
    sprintf(buffer, "TraLog%i", n + 1);
  	NTraLog[n] = gGeoManager->MakeBox(buffer,pMedSi,3.0 / 2., 1.5 / 2., thAlpide / 2.);
  	NTraLog[n]->SetVisLeaves(kTRUE);
  	NTraLog[n]->SetLineColor(kOrange);
  }
  
  TGeoCombiTrans *pMatrixd[36];
  float alpx[3] = { 1.5004, -1.5, -4.5004 };
  float offsetxz = (0.0199 + thAlpide);
  
  TGeoRotation* crot = new TGeoRotation();
  crot->RotateX(180.);
  for(Int_t n = 0; n < 18; n++){
  	float pdx;
  	if(n%3 == 0) pdx = alpx[0];
  	if(n%3 == 1) pdx = alpx[1];
  	if(n%3 == 2) pdx = alpx[2];
  	float pdy;
  	float pdz = 0.;
  	if(n/3 == 0) pdy = -1.5004 / 2.0 - 3.0 * 1.0;
  	if(n/3 == 1) pdy =  1.5004 / 2.0 - 3.0 * 1.0;
  	if(n/3 == 2) {pdy = -1.5004 / 2.0; pdx = pdx + offsetxz; pdz = pdz + offsetxz;}
  	if(n/3 == 3) {pdy = 1.5004 / 2.0; pdx = pdx + offsetxz; pdz = pdz + offsetxz;}
  	if(n/3 == 4) pdy = -1.5004 / 2.0 + 3.0 * 1.0;
  	if(n/3 == 5) pdy = 1.5004 / 2.0 + 3.0 * 1.0;
  	
  	if(n/3 == 2 || n/3 == 3) { pMatrixd[n] = new TGeoCombiTrans("", pdx, pdy, pdz, crot); }
  	else {pMatrixd[n] = new TGeoCombiTrans("", pdx, pdy, pdz, fRefRot);}
  	
  	Plane0->AddNode(NTraLog[n], 21 +n , pMatrixd[n]);
    
  }
  auto offset = 0.01595 + thAlpide / 2.0;
  auto offsetx =  -1.5;
  auto offsetz = (0. - offset - AlRingWd / 2.) ;
  auto pMatrixf00 = new TGeoCombiTrans("", offsetx, 0., -offset-offsetxz/2., fRefRot);
  offsetx = -1.5  ; 
  offsetz = (0. + offset + AlRingWd / 2.) ;
  auto pMatrixf01 = new TGeoCombiTrans("", offsetx , 0., offset + offsetxz/2., fRefRot);
  Plane0->AddNode(realring, 1, pMatrixf01);
 // realring->SetLineColor(kRed);
  Plane0->AddNode(realring, 1, pMatrixf00);
  
    TGeoRotation *alprotl = new TGeoRotation();
    alprotl->RotateX(0.);
    alprotl->RotateY(ATheta);
    alprotl->RotateZ(0.);
    

 
 
  Double_t fdx1 = 4.5*TMath::Sin(TMath::DegToRad()*Theta1_comp)+TMath::Tan(TMath::DegToRad()*alpha)*((fdz1-2.5)+ 4.5*TMath::Cos(TMath::DegToRad()*Theta1_comp));
  cout << "Pos alpide X : " << fdx1 << " cm" << endl;
  cout << "Pos alpide Z : " << fdz1 << " cm" << endl;
  cout << "Dist Alpide : " << TMath::Sqrt(fdx1*fdx1+(fdz1+shift-2.5)*(fdz1+shift-2.5)) << " cm" << endl;
  fdx1 = fdx1 + 1.5*TMath::Sin(TMath::DegToRad()*Theta1_comp);
    float fdy1 = 0.0;
  fdz1 = fdz1 - 1.5*TMath::Cos(TMath::DegToRad()*Theta1_comp) + 0.;
  
  TGeoCombiTrans *alpt00l = new TGeoCombiTrans(fdx1,fdy1,fdz1+shift,alprotl);
  
  pWorld->AddNode(Plane0, 1, GetGlobalPosition(alpt00l));

// the coordinates in the alpide, transplated from the center position which was got from the foot.
  cout << "the exit angle: " <<  TMath::ATan((fdx1 - 6*TMath::Cos(TMath::DegToRad()*Theta1))/((fdz1-2.5) +shift + 6*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;
  cout << "the exit angle (0,0): " <<  TMath::ATan((fdx1 - 6*TMath::Cos(TMath::DegToRad()*Theta1))/(fdz1 +shift + 6*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;
  cout << "the exit angle (0,5): " <<  TMath::ATan((fdx1 - 6*TMath::Cos(TMath::DegToRad()*Theta1))/(fdz1-5 +shift + 6*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;
  cout << "Theta 1 : " << Theta1 << " deg" << endl;
  cout << "Theta 2 : " << Theta2 << " deg" << endl;
  //cout << "INNER angle: " << 180 - (90  - TMath::ATan((2.5 - fdz1 - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp))/(fdx1+shift + 3*TMath::Cos(TMath::DegToRad()*Theta1)))*TMath::RadToDeg()) << endl;

   //ALPIDE angular coverage
 if (fdz1-3*TMath::Cos(TMath::DegToRad()*Theta1_comp) < 2.5){
   //cout << "extends before centre" << endl;
   cout << "the angle from z axis to back of ALPIDE: " <<  180 - TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << endl;

   cout << "angular coverage: " << alpha << " to " <<  180 - TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/(2.5 - fdz1 + 3*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " (" << 180 - alpha - TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << ") and (" << 90 - (180 -  TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg()) << ")" << endl;
 }
 else{

   //cout << "doesn't extend past centre" << endl;
      cout << "the angle from z axis to back of ALPIDE: " <<  TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;

      cout << "angular coverage: " << alpha << " to " <<  TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " ("<< TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() - alpha << ") and (" << 90 -  TMath::ATan((fdx1 + 3*TMath::Cos(TMath::DegToRad()*Theta1))/((fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << ")" << endl;
 }

 //FOOT 1 angular coverage
 if (PosZ1b-5*TMath::Cos(TMath::DegToRad()*Theta1_comp) < 2.5){
   //cout << "extends before centre" << endl;
   cout << "the angle from z axis to back of FOOT1: " <<  180 - TMath::ATan((PosX1b + 0.0075*TMath::Sin(TMath::DegToRad()*Theta1) + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp) + 0.0075*TMath::Sin(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << endl;

   cout << "angular coverage: " <<  "Theta exit : " << TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() <<" deg" << " to " <<  180 - TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(2.5 - PosZ1b + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " (" << 180 - TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() - TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << ") and (" << 90 - (180 -  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg()) << ")" << endl;
 }
 else{

   //cout << "doesn't extend past centre" << endl;
   cout << "the angle from z axis to back of FOOT1: " <<  TMath::ATan((PosX1b + 0.0075*TMath::Sin(TMath::DegToRad()*Theta1) + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift + 0.0075*TMath::Sin(TMath::DegToRad()*Theta1_comp) - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;

      cout << "angular coverage: " << TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " to " <<  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " ("<< TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() - TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << ") and (" << 90 -  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << ")" << endl;
 }
 /*
 if (PosZ1b-5*TMath::Cos(TMath::DegToRad()*Theta1_comp) < 2.5){
   //cout << "extends before centre" << endl;
   cout << "the angle from z axis to back of FOOT1: " <<  180 - TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << endl;

   cout << "angular coverage: " <<  "Theta exit : " << TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() <<" deg" << " to " <<  180 - TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(2.5 - PosZ1b + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " (" << 180 - TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() - TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg() << ") and (" << 90 - (180 -  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/(-1*((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp))))*TMath::RadToDeg()) << ")" << endl;
 }
 else{

   //cout << "doesn't extend past centre" << endl;
      cout << "the angle from z axis to back of FOOT1: " <<  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;

      cout << "angular coverage: " << TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " to " <<  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << " ("<< TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() - TMath::ATan((PosX1b - 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << ") and (" << 90 -  TMath::ATan((PosX1b + 5*TMath::Cos(TMath::DegToRad()*Theta1))/((PosZ1b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << ")" << endl;
 }
*/
 
  //FOOT 2 angular coverage
 if (PosZ2b-5*TMath::Cos(TMath::DegToRad()*Theta1_comp) < 2.5){
   //cout << "extends before centre" << endl;
   cout << "the angle from z axis to back of FOOT2: " <<  180 - TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/(-1*((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp))))*TMath::RadToDeg() << endl;

   cout << "angular coverage: " <<  "Theta exit : " << TMath::ATan((PosX2b - 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() <<" deg" << " to " <<  180 - TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/(2.5 - PosZ2b + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << " (" << 180 - TMath::ATan((PosX2b - 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() - TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/(-1*((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp))))*TMath::RadToDeg() << ") and (" << 90 - (180 -  TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/(-1*((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp))))*TMath::RadToDeg()) << ")" << endl;
 }
 else{

   //cout << "doesn't extend past centre" << endl;
      cout << "the angle from z axis to back of FOOT2: " <<  TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << endl;

      cout << "angular coverage: " << TMath::ATan((PosX2b - 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << " to " <<  TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << " ("<< TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() - TMath::ATan((PosX2b - 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5) + shift + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << ") and (" << 90 -  TMath::ATan((PosX2b + 5*TMath::Cos(TMath::DegToRad()*Theta2))/((PosZ2b-2.5)+shift - 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << ")" << endl;
 }
 cout << " ALPIDE forward distance: " << TMath::Sqrt(pow((fdx1  - 6*TMath::Sin(TMath::DegToRad()*Theta1_comp)),2) + pow((fdz1 - 2.5 + 6*TMath::Cos(TMath::DegToRad()*Theta1_comp)),2))<< endl;
 cout << " FOOT 1 forward distance: " << TMath::Sqrt(pow((PosX1b  - 5*TMath::Sin(TMath::DegToRad()*Theta1_comp)),2) + pow((PosZ1b - 2.5 + 5*TMath::Cos(TMath::DegToRad()*Theta1_comp)),2)) << endl;
 cout << " FOOT 2 forward distance: " << TMath::Sqrt(pow((PosX2b  - 5*TMath::Sin(TMath::DegToRad()*Theta2_comp)),2) + pow((PosZ2b - 2.5 + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)),2))<< endl;
 //cout << " FOOT 1 angle: " << TMath::ATan((PosX2b  - 5*TMath::Sin(TMath::DegToRad()*Theta2_comp))/(PosZ2b - 2.5 + 5*TMath::Cos(TMath::DegToRad()*Theta2_comp)))*TMath::RadToDeg() << endl;

 cout << "flank distance: " << fdx1 - 8.2*TMath::Sin(TMath::DegToRad()*Theta1_comp) - 0.644972*TMath::Cos(TMath::DegToRad()*Theta1_comp) - 0.09*TMath::Sin(TMath::DegToRad()*Theta1_comp) << endl;
 cout << "flank distance old: " << (fdx1 - 8.2*TMath::Sin(TMath::DegToRad()*Theta1_comp))*20  << endl;
 cout << "flank distance z: " << fdz1 -2.5 + 8.2*TMath::Cos(TMath::DegToRad()*Theta1_comp) << endl;
cout << "flank angle: " << TMath::ATan((fdx1 - 8.2*TMath::Sin(TMath::DegToRad()*Theta1_comp))/(fdz1 - 2.5 +8.2*TMath::Cos(TMath::DegToRad()*Theta1_comp)))*TMath::RadToDeg() << endl;
 // cout << "  ...  " << TMath::ATan((0.622)/(3.6739))*TMath::RadToDeg() << endl;
// cout << " .... " << (fdz1-2.5)+shift - 3*TMath::Cos(TMath::DegToRad()*Theta1_comp) << endl;
 cout << "ALPIDE - FOOT1 gap: " << GapAF << "cm" << endl;  
 cout << "FOOT1 - FOOT2 gap: " << GapF << "cm" << endl;

 /*-------------------------------------------------------ALPIDE RIGHT ARM--------------------------------------------------*/ 
  
  TGeoVolume* Plane1 = new TGeoVolumeAssembly("APlane1");
  TGeoRotation* brot = new TGeoRotation();
  brot->RotateZ(180.);
  TGeoRotation* drot = new TGeoRotation();
  drot->RotateX(180.);
  drot->RotateZ(180.);
  for(Int_t n = 0; n < 18; n++){
  	float pdx;
  	if(n%3 == 0) pdx = alpx[0];
  	if(n%3 == 1) pdx = alpx[1];
  	if(n%3 == 2) pdx = alpx[2];
  	float pdy;
  	float pdz = 0.;
  	if(n/3 == 0) pdy = -1.5004 / 2.0 - 3.0 * 1.0;
  	if(n/3 == 1) pdy =  1.5004 / 2.0 - 3.0 * 1.0;
  	if(n/3 == 2) {pdy = -1.5004 / 2.0; pdx = pdx - offsetxz; pdz = pdz + offsetxz;}
  	if(n/3 == 3) {pdy =  1.5004 / 2.0; pdx = pdx - offsetxz; pdz = pdz + offsetxz;}
  	if(n/3 == 4) pdy = -1.5004 / 2.0 + 3.0 * 1.0;
  	if(n/3 == 5) pdy =  1.5004 / 2.0 + 3.0 * 1.0;
  	
  	if(n/3 == 2 || n/3 == 3) { pMatrixd[n+18] = new TGeoCombiTrans("", pdx, pdy, pdz, drot); }
  	else {pMatrixd[n+18] = new TGeoCombiTrans("", pdx, pdy, pdz, brot);}
  	
  	Plane1->AddNode(NTraLog[n+18], 21 + 20 +n, pMatrixd[n+18]);
    
  }
  offset = 0.01595 + thAlpide / 2.0;
  offsetx =  1.5;
  offsetz = (0. - offset - AlRingWd / 2.) ;
  auto pMatrixf10 = new TGeoCombiTrans("", offsetx, 0., -offset-offsetxz/2., brot);
  offsetx = 1.5  ; 
  offsetz = (0. + offset + AlRingWd / 2.) ;
  auto pMatrixf11 = new TGeoCombiTrans("", offsetx , 0., offset + offsetxz/2., brot);
  Plane1->AddNode(realring, 1, pMatrixf01);
 // realring->SetLineColor(kRed);
  Plane1->AddNode(realring, 1, pMatrixf00);  
  
    TGeoRotation *alprotr = new TGeoRotation();
    alprotr->RotateX(0.);
    alprotr->RotateZ(180.);
    alprotr->RotateY(-1*ATheta);
    double fdx2 = -fdx1;
    double fdy2 = 0.0;
    double fdz2 = fdz1;
  
  TGeoCombiTrans *alpt00r = new TGeoCombiTrans(fdx2,fdy2,fdz2+shift,alprotr);
  
  pWorld->AddNode(Plane1, 1, GetGlobalPosition(alpt00r));
  
  
/*	float fdz = 0.;
	auto pMatrixd0 = new TGeoCombiTrans("", 0., 0., fdz, fRefRot);
	float x[nbsensorpermultilayer] = { 1.5004, -1.5, -4.5004, 1.5004, -1.5, -4.5004 };
	for (Int_t d = 0; d < nbsensorpermultilayer; d++)
	{
	    float y = -1.5004 / 2.0 + (d > 2 ? 1.5004 : 0.0);
	    pMatrixd0 = new TGeoCombiTrans("", x[d], y, fdz, fRefRot);
	}
  
  
  

    auto Alpide = gGeoManager->MakeBox("Alpide", pMedSi, 3.0 / 2., 1.5 / 2., thAlpide / 2.);

    sprintf(buffer, "Multilayer_%i", 1);
    TGeoVolume* Layer =
        gGeoManager->MakeBox(buffer, pMedVac, Multilayerlength / 2.0, MultilayerWd / 2.0, 0.00395 + thAlpide / 2.0);

    Alpide->SetLineColor(kOrange);
    Alpide->SetVisLeaves(kTRUE);



    TGeoVolume* Plane1 = new TGeoVolumeAssembly("APlane1");
    
    TGeoRotation* arot = new TGeoRotation();
    arot->RotateY(45.);
    TGeoRotation* brot = new TGeoRotation();
    brot->RotateZ(180.);
    brot->RotateY(-45.);

    TGeoRotation* crot = new TGeoRotation();
    crot->RotateX(180.);
    crot->RotateY(45.);
    TGeoRotation* drot = new TGeoRotation();
    drot->RotateX(180.);
    drot->RotateZ(180.);
    drot->RotateY(-45.);

    float fdx1 = 3.5;
    float fdz1 = 8.0 - 1.5;
    float offsetxz = (0.0199 + thAlpide) * cos(45. * TMath::DegToRad());
    
    
    for (Int_t n = 0; n < nbmultilayers; n++)
    {

            if (n == 0)
            {
                float fdz = 0.;
                auto pMatrixd0 = new TGeoCombiTrans("", 0., 0., fdz, fRefRot);
                float x[nbsensorpermultilayer] = { 1.5004, -1.5, -4.5004, 1.5004, -1.5, -4.5004 };
                for (Int_t d = 0; d < nbsensorpermultilayer; d++)
                {
                    float y = -1.5004 / 2.0 + (d > 2 ? 1.5004 : 0.0);
                    pMatrixd0 = new TGeoCombiTrans("", x[d], y, fdz, fRefRot);
                    Layer->AddNode(Alpide, d + 1, pMatrixd0);
                }
            }

            auto pMatrixd0 = new TGeoCombiTrans("", fdx1, -3.0 * (1.0 - n), fdz1, arot);

            if (n == 0 || n == 2)
            {
                Plane0 ->AddNode(Layer, 1 + n, pMatrixd0);
            }
            else
            {
                pMatrixd0 = new TGeoCombiTrans(
                    "", fdx1 + offsetxz , -3.0 * (1.0 - n), fdz1 + offsetxz, crot);
                Plane0 ->AddNode(Layer, 1 + n, pMatrixd0);

                auto offset = 0.01595 + thAlpide / 2.0;

                pMatrixd0 = new TGeoCombiTrans("", fdx1, 0., fdz1 + offset, arot);

                auto offsetx =  -(1.5 + offset + AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                auto offsetz = (1.5 - offset - AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                pMatrixd0 = new TGeoCombiTrans("", fdx1 + offsetx, 0., fdz1 + offsetz, arot);
                Plane0->AddNode(realring, 1, pMatrixd0);

                offsetx =  -(1.5 - offset - AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                offsetz = (1.5 + offset + AlRingWd / 2.) * cos(45 * TMath::DegToRad());

                pMatrixd0 = new TGeoCombiTrans("",  fdx1 + offsetx + offsetxz, 0., fdz1 + offsetz + offsetxz, arot);
                Plane0 ->AddNode(realring, 1, pMatrixd0);
            }
    }
        
    fdx1 = -3.5;
    fdz1 = 8.0 - 1.5; 
    for (Int_t n = 0; n < nbmultilayers; n++)
    {

            if (n == 0)
            {
                float fdz = 0.;
                auto pMatrixd0 = new TGeoCombiTrans("", 0., 0., fdz, fRefRot);
                float x[nbsensorpermultilayer] = { 1.5004, -1.5, -4.5004, 1.5004, -1.5, -4.5004 };
                for (Int_t d = 0; d < nbsensorpermultilayer; d++)
                {
                    float y = -1.5004 / 2.0 + (d > 2 ? 1.5004 : 0.0);
                    pMatrixd0 = new TGeoCombiTrans("", x[d], y, fdz, fRefRot);
                    Layer->AddNode(Alpide, d + 1, pMatrixd0);
                }
            }

            auto pMatrixd0 = new TGeoCombiTrans("", fdx1, -3.0 * (1.0 - n), fdz1, brot);

            if (n != 1)
            {
                Plane1->AddNode(Layer, 1 + n, pMatrixd0);
            }
            else
            {
                pMatrixd0 = new TGeoCombiTrans( "", fdx1 + -offsetxz, -3.0 * (1.0 - n), fdz1 + offsetxz, drot);
                Plane1->AddNode(Layer, 1 + n, pMatrixd0);
            
                auto offset = 0.01595 + thAlpide / 2.0;

                pMatrixd0 = new TGeoCombiTrans("", fdx1, 0., fdz1 + offset, brot);

                auto offsetx =  (1.5 + offset + AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                auto offsetz = (1.5 - offset - AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                pMatrixd0 = new TGeoCombiTrans("", fdx1 + offsetx, 0., fdz1 + offsetz, brot);
                Plane1->AddNode(realring, 1, pMatrixd0);

                offsetx =  (1.5 - offset - AlRingWd / 2.) * cos(45 * TMath::DegToRad());
                offsetz = (1.5 + offset + AlRingWd / 2.) * cos(45 * TMath::DegToRad());

                pMatrixd0 = new TGeoCombiTrans("", fdx1 + offsetx  -offsetxz, 0.,  fdz1 + offsetz + offsetxz, brot);
                Plane1->AddNode(realring, 1, pMatrixd0);
            }
        } 

  pWorld->AddNode(Plane0, 1, GetGlobalPosition(t00));
  pWorld->AddNode(Plane1, 1, GetGlobalPosition(t00));
*/
  
}

void create_califa_target_geo(const int index = 0)
{
    if (index == 1)
    {
        create_califa_geo_selector("s467");
    }
    else if (index == 2)
    {
        create_califa_geo_selector("s455");
    }
    else if (index == 3)
    {
        create_califa_geo_selector("s118");
    }
    else
    {
        create_califa_geo_selector();
    }
}
TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef)
{
  if (fLocalTrans == kTRUE ) {
    
    if ( ( fThetaX == 0 )  && ( fThetaY==0 )  && ( fThetaZ == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    // X axis
    Double_t xAxis[3] = { 1. , 0. , 0. };
    Double_t yAxis[3] = { 0. , 1. , 0. };
    Double_t zAxis[3] = { 0. , 0. , 1. };
    // Reference Rotation
    fRefRot = fRef->GetRotation();
    
    if (fRefRot) {
      Double_t mX[3] = {0.,0.,0.};
      Double_t mY[3] = {0.,0.,0.};
      Double_t mZ[3] = {0.,0.,0.};
      
      fRefRot->LocalToMasterVect(xAxis,mX);
      fRefRot->LocalToMasterVect(yAxis,mY);
      fRefRot->LocalToMasterVect(zAxis,mZ);
      
      Double_t a[4]={ mX[0],mX[1],mX[2], fThetaX };
      Double_t b[4]={ mY[0],mY[1],mY[2], fThetaY };
      Double_t c[4]={ mZ[0],mZ[1],mZ[2], fThetaZ };
      
      ROOT::Math::AxisAngle aX(a,a+4);
      ROOT::Math::AxisAngle aY(b,b+4);
      ROOT::Math::AxisAngle aZ(c,c+4);
      
      ROOT::Math::Rotation3D fMatX( aX );
      ROOT::Math::Rotation3D fMatY( aY );
      ROOT::Math::Rotation3D fMatZ( aZ );
      
      ROOT::Math::Rotation3D  fRotXYZ = (fMatZ * (fMatY * fMatX));
      
      //cout << fRotXYZ << endl;
      
      Double_t fRotable[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      fRotXYZ.GetComponents(
			    fRotable[0],fRotable[3],fRotable[6],
			    fRotable[1],fRotable[4],fRotable[7],
			    fRotable[2],fRotable[5],fRotable[8]
			    );
      TGeoRotation *pRot = new TGeoRotation();
      pRot->SetMatrix(fRotable);
      TGeoCombiTrans *pTmp = new TGeoCombiTrans(*fGlobalTrans,*pRot);
      
      // ne peut pas etre applique ici
      // il faut differencier trans et rot dans la multi.
      TGeoRotation rot_id;
      rot_id.SetAngles(0.0,0.0,0.0);
      
      TGeoCombiTrans c1;
      c1.SetRotation(rot_id);
      const Double_t *t = pTmp->GetTranslation();
      c1.SetTranslation(t[0],t[1],t[2]);
      
      TGeoCombiTrans c2;
      c2.SetRotation(rot_id);
      const Double_t *tt = fRefRot->GetTranslation();
      c2.SetTranslation(tt[0],tt[1],tt[2]);
      
      TGeoCombiTrans cc = c1 * c2 ;
      
      TGeoCombiTrans c3;
      c3.SetRotation(pTmp->GetRotation());
      TGeoCombiTrans c4;
      c4.SetRotation(fRefRot);
      
      TGeoCombiTrans ccc = c3 * c4;
      
      TGeoCombiTrans pGlobal;
      pGlobal.SetRotation(ccc.GetRotation());
      const Double_t *allt = cc.GetTranslation();
      pGlobal.SetTranslation(allt[0],allt[1],allt[2]);
      
      return  ( new TGeoCombiTrans( pGlobal ) );
      
    }else{
      
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      No. Ref. Transformation defined ! " << endl;
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      cannot create Local Transformation " << endl;
      return NULL;
    } //! fRefRot
    
  } else {
    // Lab Transf.
    if ( ( fPhi == 0 )  && ( fTheta==0 )  && ( fPsi == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    return ( new TGeoCombiTrans(*fGlobalTrans,*fGlobalRot) );
    
  }
}
