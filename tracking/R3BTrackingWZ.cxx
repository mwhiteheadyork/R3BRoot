/******************************************************************************
 *   Copyright (C) 2024 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2024 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// Created on 09/02/2024 by M.Xarepe

#include "R3BTrackingWZ.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BEventHeader.h"
#include "R3BLosHitData.h"
#include "R3BTofdHitData.h"
#include "R3BMwpcHitData.h"
#include "R3BFrsData.h"
#include "R3BTttxHitData.h"

#include "R3BMCTrack.h"
#include "R3BMDFWrapper.h"
#include "R3BTrack.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "R3BShared.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TCutG.h"
#include "TFile.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <R3BCoarseTimeStitch.h>
#include "THttpServer.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
R3BTrackingWZ* gMDFTrackerWZ;

R3BTrackingWZ::R3BTrackingWZ()
	: R3BTrackingWZ("TrackingWZ", 1)
{
}

R3BTrackingWZ::R3BTrackingWZ(const char* name, Int_t iVerbose)
	: FairTask(name, iVerbose)
	, fTrigger(-1)
	, fTpat(-1)
	, fNEvents(0)
	, maxevent(0)
	  , DoAlignment(false)
	, fTrackItems(new TClonesArray("R3BTrack"))
	, reference_PoQ(0.)
	, GladCurrent(-1)
	, GladReferenceCurrent(-1)
	, FiberTimeMin(-1)
	, FiberTimeMax(-1)
	, FiberEnergyMin(-1)
	, FiberEnergyMax(-1)
	  , fHeader(nullptr)
{
}

R3BTrackingWZ::~R3BTrackingWZ()
{
	if (fTrackItems){
		delete fTrackItems;}
}

InitStatus R3BTrackingWZ::Init()
{
	LOG(info) << "R3BTrackingWZ::Init()";
	FairRootManager* mgr = FairRootManager::Instance();
	if (NULL == mgr)
	{
		LOG(fatal) << "FairRootManager not found";
	}
	fTimeStitch = new R3BCoarseTimeStitch();
	fTimeStitch->SetClockTDC150();
	fHeader = (R3BEventHeader*)mgr->GetObject("EventHeader.");
	//fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
	if (!fHeader)
	{
		LOG(warn) << "R3BTrackingWZ::Init() EventHeader. not found";
	}
	// Reading all detector branches
	cout << "\nDET_MAX = " << DET_MAX << endl;
	assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
	LOG(info) << "Reading " << NOF_FIB_DET << " fiber detectors";
	for (int det = 0; det < DET_MAX; det++)
	{
		fDataItems.push_back((TClonesArray*)mgr->GetObject(Form("%s", fDetectorNames[det])));
		if (NULL == fDataItems.at(det))
		{
			R3BLOG(fatal, Form("\n\n Cannot find tree branch %s \n\n", fDetectorNames[det]));
		}
	}
	// check if all cuts are properly set
	if (GladCurrent < 0 || GladReferenceCurrent < 0 || 
			FiberEnergyMin < 0 || FiberEnergyMax < 0)
	{
		R3BLOG(fatal, Form(" Some cuts are not set or negative values are used\n\n"));
	}
	// Initializing all MDF functions
	LOG(info) << "Reading MDF function for TX0";
	MDF_TX0 = new R3BMDFWrapper(MDF_TX0_filename.Data());

	LOG(info) << "Reading MDF function for FlightPath";
	MDF_FlightPath = new R3BMDFWrapper(MDF_FlightPath_filename.Data());

	LOG(info) << "Reading MDF function for TY0";
	MDF_TY0 = new R3BMDFWrapper(MDF_TY0_filename.Data());

	LOG(info) << "Reading MDF function for TX1";
	MDF_TX1 = new R3BMDFWrapper(MDF_TX1_filename.Data());

	LOG(info) << "Reading MDF function for TY1";
	MDF_TY1 = new R3BMDFWrapper(MDF_TY1_filename.Data());

	LOG(info) << "Reading MDF function for PoQ";
	MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_filename.Data());


	//Read output from the vertex macro
	// linking to global pointer (needed by alignment)
	gMDFTrackerWZ = this;
	// Request storage of R3BTrack data in the output tree
	mgr->Register("MDFTracks", "MDFTracks data", fTrackItems, kTRUE);

	//online

    	TFolder* mainfol = new TFolder("tracker", "tracker_info");
	FairRunOnline* run = FairRunOnline::Instance();
//	run->GetHttpServer()->Register("", this);
	run->AddObject(mainfol);

//	run->GetHttpServer()->RegisterCommand("Reset_tracker", Form("/Objects/%s/->Reset_Tracker_Histo()", GetName()));

	trackerCanvas = new TCanvas("tracker_Canvas", "trackerCanvas");
	AoQ_Vs_Q_TOFD =
		R3B::root_owned<TH2F>("AoQ_Vs_Q_tofd", "AoQ_Vs_Q_tofd", 1000, 1., 5, 1200, 0, 12);

	AoQ_Vs_Q_TOFD->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD->GetYaxis()->SetTitle("tofd Q");
	trackerCanvas->cd();
	AoQ_Vs_Q_TOFD->Draw("COLZ");
	mainfol->Add(trackerCanvas);
	
	AoQ_Vs_Q_TOFD_fi31 =
		R3B::root_owned<TH2F>("AoQ_Vs_Q_tofd_fi31", "AoQ_Vs_Q_tofd_fi31", 1000, 1., 5, 1200, 0, 12);

	AoQ_Vs_Q_TOFD_fi31->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD_fi31->GetYaxis()->SetTitle("tofd Q");
	trackerCanvas->cd();
	AoQ_Vs_Q_TOFD_fi31->Draw("COLZ");

	AoQ_Vs_Q_TOFD_fi33 =
		R3B::root_owned<TH2F>("AoQ_Vs_Q_tofd_fi33", "AoQ_Vs_Q_tofd_fi33", 1000, 1., 5, 1200, 0, 12);

	AoQ_Vs_Q_TOFD_fi33->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD_fi33->GetYaxis()->SetTitle("tofd Q");
	trackerCanvas->cd();
	AoQ_Vs_Q_TOFD_fi33->Draw("COLZ");

	AoQ_Vs_Q_TOFD_bothfi =
		R3B::root_owned<TH2F>("AoQ_Vs_Q_tofd_bothfi", "AoQ_Vs_Q_tofd_bothfi", 1000, 1., 5, 1200, 0, 12);

	AoQ_Vs_Q_TOFD_bothfi->GetXaxis()->SetTitle("AoQ");
	AoQ_Vs_Q_TOFD_bothfi->GetYaxis()->SetTitle("tofd Q");
	trackerCanvas->cd();
	AoQ_Vs_Q_TOFD_bothfi->Draw("COLZ");

	Track_mult_bothfi = R3B::root_owned<TH1F>("bothfi","bothfi", 10, 0, 10);
	trackerCanvas->cd();
	Track_mult_bothfi->Draw("COLZ");
	Track_mult_fi31 = R3B::root_owned<TH1F>("fi31","fi31", 10, 0, 10);
	trackerCanvas->cd();
	Track_mult_fi31->Draw("COLZ");
	Track_mult_fi33 = R3B::root_owned<TH1F>("fi33","fi33", 10, 0, 10);
	trackerCanvas->cd();
	Track_mult_fi33->Draw("COLZ");

	return kSUCCESS; 
}

void R3BTrackingWZ::Reset_Tracker_Histo(){
	AoQ_Vs_Q_TOFD->Reset();
	return;
}

void R3BTrackingWZ::Exec(Option_t* option)
{
	if (fNEvents / 1000. == (int)fNEvents / 1000)
		std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
			<< " %) " << std::flush;
	FairRootManager* mgr = FairRootManager::Instance();
	R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");
	fNEvents += 1;
	is_good_event = false;
	Tpat = fHeader->GetTpat();//vairable in the output tree
				

	if((Tpat & 0xf000)){
		return;
	}

	if (fTpat != -1){
		if ((Tpat & fTpat)!=fTpat){
			return;
		}
	}

	mul_los=-999;
	mul_m0=-999;
	mul_m1=-999;
	mul_f32=-999;
	mul_f30=-999;
	mul_f31=-999;
	mul_f33=-999;
	mul_tofd=-999;
	cond=false;

	fi31_only = false;
	fi33_only = false;
	bothfi = false;
	fi31_track = false;
	fi33_track = false;

	mul_los   = fDataItems[LOS_DATA]->GetEntriesFast();
	mul_m0   = fDataItems[MWPC0_HITDATA]->GetEntriesFast();
	mul_m1   = fDataItems[MWPC1_HITDATA]->GetEntriesFast();
	//mul_tttx = fDataItems[TTTX_HITDATA]->GetEntriesFast();
	mul_f32  = fDataItems[DET_FI32]->GetEntriesFast();
	mul_f30  = fDataItems[DET_FI30]->GetEntriesFast();
	mul_f31  = fDataItems[DET_FI31]->GetEntriesFast();
	mul_f33  = fDataItems[DET_FI33]->GetEntriesFast();
	mul_tofd = fDataItems[DET_TOFD]->GetEntriesFast();
	//if(mul_los!=1) return;
	a++;
	//if(mul_tttx<1) return;
	if(mul_tofd<1) return;
	b++;
	//if(mul_m0!=1)return;
	//if(mul_m1!=1) return;
	c++;
	if(mul_f32<1 || mul_f30<1 || (mul_f31==0 && mul_f33==0)) return;
	d++;
	//if(mul_foot<1) return;//for now take only mul=1 in mwpcs
	//FRS data
	//auto frs_DataItems = fDataItems.at(FRS_DATA);
	//if(frs_DataItems->GetEntriesFast() < 1) return; 
	//auto frs_data = (R3BFrsData*)frs_DataItems->At(0);
	//if(frs_data->GetBrho()<17 || frs_data->GetBrho()>18) return;
	//if(frs_data->GetAq()<2.675 || frs_data->GetAq()>2.694) return;
	//if(frs_data->GetZ()<5.2 || frs_data->GetZ()>6.7) return;
	//cout << "\nGood event!\n";

	//------ Get TOFD data 
	R3BTofdHitData* tofd_hit{};
	int mul_tofd1=0;
	double tofdq_temp =0;
	double tttxq_temp =0;
	bool is_good_tofd = false;
	bool is_good_tttx = false;
	for (auto i = 0; i < fDataItems[DET_TOFD]->GetEntriesFast(); ++i)
	{
		tofd_hit = static_cast<R3BTofdHitData*>(fDataItems[DET_TOFD]->At(i));
		if (tofd_hit->GetDetId() == 1) // only hits from first plane, add Z later
		{
			//if (tofd_hit->GetTime()<0 || tofd_hit->GetTime()>-200){
				is_good_tofd = true;
				mul_tofd1++;
				tofdq_temp=tofd_hit->GetEloss();
			//}
			//break;
		}
	}
//	for (int ihit = 0; ihit < fDataItems[TTTX_HITDATA]->GetEntriesFast(); ihit++)
//	{
//		auto hit = dynamic_cast<R3BTttxHitData*>(fDataItems[TTTX_HITDATA]->At(ihit));
//		auto idet = hit->GetDetID();
//		if(idet==1){
//			tttxq_temp = hit->GetChargeZ();
//			is_good_tttx = true;
//		}
//	}
//	if(!is_good_tttx){ 
//		return;
//	}
	if(!is_good_tofd){ 
		return;
	}
	//if(mul_tofd1!=1){ 
	//	return;
	//}
	e++;
	if(!MakeIncomingTracks()){ 
		return;//at least one good track candidate in FOOT
	}
	if (tracks_in.size()!=1){
		cout << " BAD MWPC Track? " << endl;
		return;
	}
	f++;
	if(!MakeOutgoingTracks()){ 
		return;//at least one good track candidate in Fibers
	}

	if (bothfi){
		if (!(tracks_out.size()==1 || (tracks_out.size()==2 && fi31_track && fi33_track))){
			return;
		}
	}
	if ((fi31_only || fi33_only)){
		if (tracks_out.size()!=1){
			return;
		}
	}
			
	g++;
	is_good_event = true;
	cond=true;
	int counter = 0;
	double delta_TX1, delta_TX0, delta_TY0;
	for (auto & tin : tracks_in){
		for (auto & tout : tracks_out){
			counter++;
			//preserve the order, it is expected by the MDF function!
			mdf_data[0] = tin.mw1_x;
			mdf_data[1] = tin.mw1_y;
			mdf_data[2] = tin.mw1_z;
			//mdf_data[0] = 0.;
			//mdf_data[1] = 0.;
			//mdf_data[2] = tin.mw1_z;
		//	mdf_data[3] = 0.0;/*(tin.mw0_x - tin.mw1_x)/(tin.mw0_z - tin.mw1_z);*/
		//	mdf_data[4] = 0;/*(tin.mw0_y - tin.mw1_y)/(tin.mw0_z - tin.mw1_z);*/
			mdf_data[3] = tout.f32_x;
			mdf_data[4] = tout.f32_z;
			mdf_data[5] = (tout.last_x - tout.f32_x)/(tout.last_z - tout.f32_z);
			mdf_data[6] = (tout.f30_y  - tin.mw0_y)/(tout.f30_z - tin.mw0_z);
			//mdf_data[6] = (tout.f30_y  - tin.mw1_y)/(tout.f30_z - tin.mw1_z);
			// Calculate all required MDF values

			flight_p = MDF_FlightPath->MDF(mdf_data);
			poq = MDF_PoQ->MDF(mdf_data) * GladCurrent / GladReferenceCurrent;
			tx0 = MDF_TX0->MDF(mdf_data);
			tx1 = MDF_TX1->MDF(mdf_data);
			ty0 = MDF_TY0->MDF(mdf_data);
			ty1 = MDF_TY1->MDF(mdf_data);
			tof = flight_p / FRS_BETA / SPEED_OF_LIGHT;
		//	hh = poq;
		//	ii = tx0;
		//	jj = ty0;

			beta = flight_p / tof / SPEED_OF_LIGHT;
			gamma = 1. / sqrt(1 - pow(beta, 2));
			maoz = poq / beta / gamma / AMU;
		//	kk = maoz;
			TVector3 vec_PoQ(tx0, ty0, 1);
			vec_PoQ.SetMag(poq);

			if(counter == 1){
				AoQ_Vs_Q_TOFD->Fill(maoz,tofdq_temp);

				if (fi31_only){
					AoQ_Vs_Q_TOFD_fi31->Fill(maoz,tofdq_temp);
					fi31_counts++;
				}
				if (fi33_only){
					AoQ_Vs_Q_TOFD_fi33->Fill(maoz,tofdq_temp);
					fi33_counts++;
				}
				if (bothfi){
					AoQ_Vs_Q_TOFD_bothfi->Fill(maoz,tofdq_temp);
					bothfi_counts++;
				}
			}

			AddTrackData(tin.mw1_x, tin.mw1_y, tin.mw1_z, vec_PoQ, tofdq_temp, maoz); // chix, chiy, quality
			if (fi31_only){
				Track_mult_fi31->Fill(counter);
			}
			if (fi33_only){
				Track_mult_fi33->Fill(counter);
			}
			if (bothfi){
				Track_mult_bothfi->Fill(counter);
			}
		}
	}
	return;
}

void R3BTrackingWZ::FinishEvent()
{
	for (auto& DataItem : fDataItems)
	{
		DataItem->Clear();
	}
	fTrackItems->Clear();
	if (fNEvents / 1000. == (int)fNEvents / 1000)
	cout << "finish event " << a << " tofd > 1: " << b << "mwpc1 hit: " << c << " 4x fiber hit: "  << d << " good tofd: " << e << " 1 incoming track: " << f << " 1 outgoing track: " << g <<" "<< hh <<" "<< ii  <<" "<< jj << " " << kk <<endl;
	//cout << "FI31 events: " << fi31_counts << endl;
	//cout << "FI33 events: " << fi33_counts << endl;
	//cout << "Both FI events: " << bothfi_counts << endl;
	
}

void R3BTrackingWZ::FinishTask()
{
	LOG(info) << "Processed " << fNEvents << " events\n\n";
	AoQ_Vs_Q_TOFD->Write();
	AoQ_Vs_Q_TOFD_fi33->Write();
	AoQ_Vs_Q_TOFD_fi31->Write();
	AoQ_Vs_Q_TOFD_bothfi->Write();
	Track_mult_fi31->Write();
	Track_mult_fi33->Write();
	Track_mult_bothfi->Write();
	cout << "ending " << a << " " << b << " " << c << " "  << d << " " << e << " " << f << " " << g << endl;

	//cout<<"WRITE"<<endl;
}
/*

bool R3BTrackingWZ::MakeIncomingTracks()
{
	tracks_in.clear();
	TVector3 vertex_mwpc;
	Track tr;
	//Get MWPC hits, for now only first hit
	auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
	m1_point.SetXYZ(m1_hit->GetX()*0.1, m1_hit->GetY()*0.1, 0.);//cm

	TransformPoint(m1_point, &m1_angles, &m1_position);//lab
	tr.mw1_x   = m1_point.X();
	tr.mw1_y   = m1_point.Y();
	tr.mw1_z   = m1_point.Z();
	if(isnan(m1_point.X()) || isnan(m1_point.Y()) || isnan(m1_point.Z())){
		cout << "NAN MWPC TRACK? " << endl;
		return false;
	}
	tracks_in.push_back(tr);
	return true;
}
*/
bool R3BTrackingWZ::MakeIncomingTracks()
{
	tracks_in.clear();
	TVector3 vertex_mwpc;
	Track tr;
	//Get MWPC hits, for now only first hit
	//auto m1_hit = static_cast<R3BMwpcHitData*>(fDataItems[MWPC1_HITDATA]->At(0));
	m1_point.SetXYZ(-0.44, -0.6, 0.);//cm

	TransformPoint(m1_point, &m1_angles, &m1_position);//lab
	tr.mw1_x   = m1_point.X();
	tr.mw1_y   = m1_point.Y();
	tr.mw1_z   = m1_point.Z();
	if(isnan(m1_point.X()) || isnan(m1_point.Y()) || isnan(m1_point.Z())){
		cout << "NAN MWPC TRACK? " << endl;
		return false;
	}
	tracks_in.push_back(tr);
	return true;
}

bool R3BTrackingWZ::IsGoodFiberHit(R3BFiberMAPMTHitData* fhit)
{
	if((fhit->GetEloss() > FiberEnergyMin) && (fhit->GetEloss() < FiberEnergyMax) && 
			(fhit->GetTime() < 20000 && fhit->GetTime()>(-20000) )
	  )
		return true;
	else 
		return false;
}
bool R3BTrackingWZ::MakeOutgoingTracks()
{
	if(fDataItems[DET_FI32]->GetEntriesFast() == 0 || fDataItems[DET_FI30]->GetEntriesFast() == 0 || 
			(fDataItems[DET_FI33]->GetEntriesFast() == 0 && fDataItems[DET_FI31]->GetEntriesFast()==0) ) 
		return false;

	
	tracks_out.clear();
	Track tr;
	double angle_out, f30_slope, f30_offset, track_slope, track_offset;
	std::vector<double> f32x;
	std::vector<double> f30y;
	std::vector<double> flastx;
	std::vector<double> flast2x;
	std::vector<double> f32e;
	std::vector<double> f30e;
	std::vector<double> flaste;
	std::vector<double> flast2e;
	TVector3 f30_edge[2];//to extract z and x in f30
   	// fRangeVftx(8192 * 1000. / 200.)        // ns
    	// fRangeClockTDC(4096 * 1000. / 250.)    // ns
    	// fRangeClockTDC150(4096 * 1000. / 150.) // ns
	//
	if (fDataItems[DET_FI33]->GetEntriesFast()==0 && fDataItems[DET_FI31]->GetEntriesFast()>0)
		fi31_only = true;
	if (fDataItems[DET_FI33]->GetEntriesFast()>0 && fDataItems[DET_FI31]->GetEntriesFast()==0)
		fi33_only = true;
	if (fDataItems[DET_FI33]->GetEntriesFast()>0 && fDataItems[DET_FI31]->GetEntriesFast()>0)
		bothfi = true;
	
	//if (fDataItems[DET_FI33]->GetEntriesFast()==0 || fDataItems[DET_FI31]->GetEntriesFast()>0)
	//	return 0;


	for (auto i=0; i<fDataItems[DET_FI32]->GetEntriesFast(); ++i)
	{
		auto f32 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI32]->At(i));
		if(!IsGoodFiberHit(f32)) continue;
	//	double fitime = fTimeStitch->GetTime(f32->GetTime_ns() - fHeader->GetTStart(), "clocktdc150", "vftx");
		double fitime =f32->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{	f32x.push_back(f32->GetX());
			f32e.push_back(f32->GetEloss());
		}
	
	}
	for (auto i=0; i<fDataItems[DET_FI30]->GetEntriesFast(); ++i)
	{
		auto f30 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI30]->At(i));
		if(!IsGoodFiberHit(f30)) continue;
	//	double fitime = fTimeStitch->GetTime(f30->GetTime_ns() - fHeader->GetTStart(), "clocktdc150", "vftx");
		double fitime =f30->GetTime_ns() - fHeader->GetTStart();
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			f30y.push_back(f30->GetY());
			f30e.push_back(f30->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI33]->GetEntriesFast(); ++i)
	{
		auto f33 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI33]->At(i));
		if(!IsGoodFiberHit(f33)) continue;

	//	double fitime = fTimeStitch->GetTime(f33->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		double fitime =f33->GetTime_ns() - fHeader->GetTStart();
	//	jj = fitime;
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			flastx.push_back(f33->GetX());
			flaste.push_back(f33->GetEloss());
		}
	}
	for (auto i=0; i<fDataItems[DET_FI31]->GetEntriesFast(); ++i)
	{
		auto f31 = static_cast<R3BFiberMAPMTHitData*>(fDataItems[DET_FI31]->At(i));
		if(!IsGoodFiberHit(f31)) continue;
	//	double fitime = fTimeStitch->GetTime(f31->GetTime_ns() - fHeader->GetTStart(), "clocktdc", "vftx");
		double fitime =f31->GetTime_ns() - fHeader->GetTStart();
	//	kk = fitime;
		if((fitime > FiberTimeMin && fitime < FiberTimeMax))
		{
			flast2x.push_back(f31->GetX());
			flast2e.push_back(f31->GetEloss());
		}
	}
	
	if(f32x.size() < 1 || f30y.size() < 1 || (flastx.size() < 1 && flast2x.size() < 1))
		return false;
	Double_t cluster2[f32x.size()][f32x.size()][2];
	Bool_t set2[f32x.size()];
	Int_t num_hit2[f32x.size()];
	Int_t hit_clust_id2[f32x.size()];
	Int_t num_clust2 = 0;
	for(Int_t i = 0; i < f32x.size(); i++)
	{
		set2[i] = false;
		num_hit2[i] = 0;
		hit_clust_id2[i] = 0;
		for(Int_t j = 0; j < f32x.size(); j++)
		{
			cluster2[i][j][0] = 0.;
			cluster2[i][j][1] = 0.;
		}
	}

	for(Int_t i = 0; i < f32x.size(); i++)
	{
		if(!set2[i])
		{
			cluster2[num_clust2][num_hit2[num_clust2]][0] = f32x[i];
			cluster2[num_clust2][num_hit2[num_clust2]][1] = f32e[i];
			num_hit2[num_clust2]++;
			set2[i] = true;
			hit_clust_id2[i] = num_clust2;
			num_clust2++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit2[hit_clust_id2[i]]; k++)
		{
			for(Int_t j = 0; j < f32x.size(); j++)
			{
				if(set2[j])
					continue;

				if(fabs(f32x[j]-cluster2[hit_clust_id2[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id2[i];
					cluster2[id][num_hit2[id]][0] = f32x[j];
					cluster2[id][num_hit2[id]][1] = f32e[j];
					num_hit2[id]++;
					set2[j] = true;
				}
			}
		}
	}
	Double_t cluster0[f30y.size()][f30y.size()][2];
	Bool_t set0[f30y.size()];
	Int_t num_hit0[f30y.size()];
	Int_t hit_clust_id0[f30y.size()];
	Int_t num_clust0 = 0;
	for(Int_t i = 0; i < f30y.size(); i++)
	{
		set0[i] = false;
		num_hit0[i] = 0;
		hit_clust_id0[i] = 0;
		for(Int_t j = 0; j < f30y.size(); j++)
		{
			cluster0[i][j][0] = 0.;
			cluster0[i][j][1] = 0.;
		}
	}

	for(Int_t i = 0; i < f30y.size(); i++)
	{
		if(!set0[i])
		{
			cluster0[num_clust0][num_hit0[num_clust0]][0] = f30y[i];
			cluster0[num_clust0][num_hit0[num_clust0]][1] = f30e[i];
			num_hit0[num_clust0]++;
			set0[i] = true;
			hit_clust_id0[i] = num_clust0;
			num_clust0++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit0[hit_clust_id0[i]]; k++)
		{
			for(Int_t j = 0; j < f30y.size(); j++)
			{
				if(set0[j])
					continue;

				if(fabs(f30y[j]-cluster0[hit_clust_id0[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id0[i];
					cluster0[id][num_hit0[id]][0] = f30y[j];
					cluster0[id][num_hit0[id]][1] = f30e[j];
					num_hit0[id]++;
					set0[j] = true;
				}
			}
		}
	}
	Double_t cluster3[flastx.size()][flastx.size()][2];
	Bool_t set3[flastx.size()];
	Int_t num_hit3[flastx.size()];
	Int_t hit_clust_id3[flastx.size()];
	Int_t num_clust3 = 0;
	for(Int_t i = 0; i < flastx.size(); i++)
	{
		set3[i] = false;
		num_hit3[i] = 0;
		hit_clust_id3[i] = 0;
		for(Int_t j = 0; j < flastx.size(); j++)
		{
			cluster3[i][j][0] = 0.;
			cluster3[i][j][1] = 0.;
		}
	}

	for(Int_t i = 0; i < flastx.size(); i++)
	{
		if(!set3[i])
		{
			cluster3[num_clust3][num_hit3[num_clust3]][0] = flastx[i];
			cluster3[num_clust3][num_hit3[num_clust3]][1] = flaste[i];
			num_hit3[num_clust3]++;
			set3[i] = true;
			hit_clust_id3[i] = num_clust3;
			num_clust3++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit3[hit_clust_id3[i]]; k++)
		{
			for(Int_t j = 0; j < flastx.size(); j++)
			{
				if(set3[j])
					continue;

				if(fabs(flastx[j]-cluster3[hit_clust_id3[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id3[i];
					cluster3[id][num_hit3[id]][0] = flastx[j];
					cluster3[id][num_hit3[id]][1] = flaste[j];
					num_hit3[id]++;
					set3[j] = true;
				}
			}
		}
	}
	Double_t cluster1[flast2x.size()][flast2x.size()][2];
	Bool_t set1[flast2x.size()];
	Int_t num_hit1[flast2x.size()];
	Int_t hit_clust_id1[flast2x.size()];
	Int_t num_clust1 = 0;
	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		set1[i] = false;
		num_hit1[i] = 0;
		hit_clust_id1[i] = 0;
		for(Int_t j = 0; j < flast2x.size(); j++)
		{
			cluster1[i][j][0] = 0.;
			cluster1[i][j][1] = 0.;
		}
	}

	for(Int_t i = 0; i < flast2x.size(); i++)
	{
		if(!set1[i])
		{
			cluster1[num_clust1][num_hit1[num_clust1]][0] = flast2x[i];
			cluster1[num_clust1][num_hit1[num_clust1]][1] = flast2e[i];
			num_hit1[num_clust1]++;
			set1[i] = true;
			hit_clust_id1[i] = num_clust1;
			num_clust1++;
		}
		else
			continue;
		for(Int_t k = 0; k < num_hit1[hit_clust_id1[i]]; k++)
		{
			for(Int_t j = 0; j < flast2x.size(); j++)
			{
				if(set1[j])
					continue;

				if(fabs(flast2x[j]-cluster1[hit_clust_id1[i]][k][0]) < .25)
				{
					Int_t id = hit_clust_id1[i];
					cluster1[id][num_hit1[id]][0] = flast2x[j];
					cluster1[id][num_hit1[id]][1] = flast2e[j];
					num_hit1[id]++;
					set1[j] = true;
				}
			}
		}
	}
	Double_t fi32xcl[num_clust2];
	Double_t fi30ycl[num_clust0];
	Double_t filastxcl[num_clust3];
	Double_t filast2xcl[num_clust1];

	for(Int_t i = 0; i < num_clust2; i++)
	{
		fi32xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit2[i]; j++)
		{
			fi32xcl[i] += cluster2[i][j][0]*cluster2[i][j][1];
			sum_ener +=cluster2[i][j][1];
		}
		if(sum_ener > 0.)
			fi32xcl[i] = (double)fi32xcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust0; i++)
	{
		fi30ycl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit0[i]; j++)
		{
			fi30ycl[i] += cluster0[i][j][0]*cluster0[i][j][1];
			sum_ener +=cluster0[i][j][1];
		}
		if(sum_ener > 0.)
			fi30ycl[i] = (double)fi30ycl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust3; i++)
	{
		filastxcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit3[i]; j++)
		{
			filastxcl[i] += cluster3[i][j][0]*cluster3[i][j][1];
			sum_ener +=cluster3[i][j][1];
		}
		if(sum_ener > 0.)
			filastxcl[i] = (double)filastxcl[i]/sum_ener;
	}
	for(Int_t i = 0; i < num_clust1; i++)
	{
		filast2xcl[i] = 0.;
		Double_t sum_ener = 0.;
		for(Int_t j = 0; j < num_hit1[i]; j++)
		{
			filast2xcl[i] += cluster1[i][j][0]*cluster1[i][j][1];
			sum_ener +=cluster1[i][j][1];
		}
		if(sum_ener > 0.)
			filast2xcl[i] = (double)filast2xcl[i]/sum_ener;
	}
	for(auto i = 0; i < num_clust2; i ++)
	{
		f32_point.SetXYZ(fi32xcl[i], 0, 0); //cm
		TransformPoint(f32_point, &f32_angles, &f32_position);
		tr.f32_x = f32_point.X();
		tr.f32_z = f32_point.Z();
		for (auto j=0; j<num_clust0; ++j)
		{
			f30_point.SetXYZ(0,fi30ycl[j],0); //cm
			TransformPoint(f30_point, &f30_angles, &f30_position);
			tr.f30_y = f30_point.Y();
			//make combination with every hit in fibers 33 and 31
			for (auto k = 0; k<num_clust3; ++k)
			{
				flast_point.SetXYZ(filastxcl[k], 0, 0); //cm
				TransformPoint(flast_point, &f33_angles, &f33_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &f30_angles, &f30_position);
				TransformPoint(f30_edge[1], &f30_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				tracks_out.push_back(tr);
				fi31_track = true;
			}

			for (auto k = 0; k<num_clust1; ++k)
			{
				flast_point.SetXYZ(filast2xcl[k], 0, 0); //cm
				TransformPoint(flast_point, &f31_angles, &f31_position);
				tr.last_x = flast_point.X();
				tr.last_z = flast_point.Z();
				angle_out = TMath::ATan((tr.last_x - tr.f32_x)/(tr.last_z - tr.f32_z)) * TMath::RadToDeg();
				//if(angle_out>(-10.) || angle_out<(-18.)) continue;
				// We need to extrapolate Z position in f30 because it was used for Y measurement
				// Define two (X,Z) points on the f30 plane:
				//Now track every combination of upstream and downstream tracks 
				f30_edge[0].SetXYZ(-1, 0, 0);
				f30_edge[1].SetXYZ(1, 0, 0);
				TransformPoint(f30_edge[0], &f30_angles, &f30_position);
				TransformPoint(f30_edge[1], &f30_angles, &f30_position);
				// Parameterize f30 plane
				f30_slope = (f30_edge[1].X() - f30_edge[0].X()) / (f30_edge[1].Z() - f30_edge[0].Z());
				f30_offset = f30_edge[0].X() - f30_slope * f30_edge[0].Z();
				track_slope  = (tr.last_x - tr.f32_x) / (tr.last_z - tr.f32_z);
				track_offset = (tr.last_x - track_slope * tr.last_z);
				// Extrapolate final X and Z position in f30
				tr.f30_z = ((track_offset - f30_offset) / (f30_slope - track_slope));// extrapolated
				tr.f30_x = (track_slope * tr.f30_z + track_offset);// extrapolated
				tracks_out.push_back(tr);
				fi33_track = true;
			}
		}
	}
	if(tracks_out.empty()) return false;
	return true;
}

void R3BTrackingWZ::TransformPoint(TVector3& point, TVector3* rot, TVector3* trans)
{
	r.SetToIdentity();
	// First Euler rotation around Y axis
	r.RotateY(rot->Y());
	// get local X axis after first rotation
	v3_localX.SetMagThetaPhi(1, r.ThetaX(), r.PhiX());
	// Second Euler rotation around local X axis
	r.Rotate(rot->X(), v3_localX);
	// get local Z axis after second rotation
	v3_localZ.SetMagThetaPhi(1, r.ThetaZ(), r.PhiZ());
	// final rotation around local Z axis
	r.Rotate(rot->Z(), v3_localZ);
	point.Transform(r);
	point += (*trans);
	return;
}

R3BTrack* R3BTrackingWZ::AddTrackData(double x, double y, double z, TVector3 poq_vec, Double_t charge, Double_t aoz)
{
	// Filling output track info
	TClonesArray& clref = *fTrackItems;
	Int_t size = clref.GetEntriesFast();
	return new (clref[size]) R3BTrack(x, y, z, poq_vec.X(),  poq_vec.Y(), poq_vec.Z(), charge, aoz, 0., 0., 0);
}

ClassImp(R3BTrackingWZ);
