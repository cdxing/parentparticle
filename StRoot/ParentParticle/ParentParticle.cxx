#include "ParentParticle.h"
#include "StRoot/CutManager/CutManager.h"
#include "StRoot/ConstManager/ConstManager.h"
//#include "../HistManager/HistManager.h"
#include "StRoot/FlowCorrection/FlowCorrection.h"
#include "StRoot/HistManager/HistManager.h" // Load STARLibrary header files
#include "StRoot/EpProManager/EpProManager.h" // Load STARLibrary header files
#include "StRoot/FlowV0/FlowV0.h"
#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoBbcHit.h"
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/Run/run.h"
#include "StThreeVectorF.hh"
//IClass header files
#include "StRoot/IClasses/IEventPlane.h"
#include "StRoot/IClasses/IEvent.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
ClassImp(ParentParticle)

    StRefMultCorr* ParentParticle::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    ParentParticle::ParentParticle(const char* name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName)
: StMaker(name)
{
    configs.read(configFileName);
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;

    mOutPut_Phi=Form("%s_%s.root",jobid,ConstManager::MixEvent[configs.Flag_ME].Data());
}

//----------------------------------------------------------------------------- 
ParentParticle::~ParentParticle()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t ParentParticle::Init() 
{
    if(!mRefMultCorr)
    {
        mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    }
    // qapid
    mCutManager = new CutManager(configs);
    //mHistManager = new HistManager();
    //mHistManager->InitQAPID();

    // eventplane 
    // EPD
    mEpdGeom = new StEpdGeom();
    mEpProManager = new EpProManager();
    mEpProManager->InitEP();
    //const int numSubEvents = 3;
    //char subEventModes[_numSubEvents] = {'r','e','e'};

    // parentparticle
    // phi-meson
    mFlowV0 = new FlowV0(configs);
    mFlowCorrection = new FlowCorrection(configs);
    mFlowCorrection->InitReCenterCorrection(configs);
    mFlowCorrection->InitShiftCorrection(configs);


    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mFile_Phi->cd();
    mFlowV0->InitPhi();
    //TFile *mRecenteringInputFile = TFile::Open("/star/u/dchen/ana/19gev_2019/eventplane/test_recenterpar.root"); // test
    //TFile *mRecenteringInputFile = TFile::Open("/star/u/dchen/ana/19gev_2019/recenter/test_recenterpar.root"); // test
    //TFile *mRecenteringInputFile = TFile::Open("/star/u/dchen/ana/19gev_2019/shift/test_recenterpar_shift.root"); // more test
    TFile *mRecenteringInputFile = TFile::Open("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_recenterpar.root"); // pwg
    if (mRecenteringInputFile->IsZombie()) {
    	std::cout << "Error opening file with Q vector /ecentering histograms" << std::endl;
    	std::cout << "No recentering in this analysis." << std::endl;
        for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
        {
            getep_sub_recen[0][sub]     = 0;
            getep_sub_recen[1][sub]     = 0;
            
	    getep_sub_wt_recen[0][sub]     = 0;
            getep_sub_wt_recen[1][sub]     = 0;

	    getep_sub_tpc_recen[0][sub]     = 0;
            getep_sub_tpc_recen[1][sub]     = 0;
        } // event plane Psi histograms
        getep_full_recen[0]     = 0;
        getep_full_recen[1]     = 0;
       // sub w/ eta weight 
	getep_full_wt_recen[0]     = 0;
        getep_full_wt_recen[1]     = 0;
       // tpc 
	getep_full_tpc_recen[0]     = 0;
        getep_full_tpc_recen[1]     = 0;
    }
    else{
      for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
      {
          getep_sub_recen[0][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1x_epd_ABCD_sub_%d",sub + 1));
          getep_sub_recen[1][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1y_epd_ABCD_sub_%d",sub + 1));

       // sub w/ eta weight 
          getep_sub_wt_recen[0][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1x_epd_ABCD_wt_sub_%d",sub + 1));
          getep_sub_wt_recen[1][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq1y_epd_ABCD_wt_sub_%d",sub + 1));

       //  tpc  
          getep_sub_tpc_recen[0][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq2x_tpc_AB_sub_%d",sub + 1));
          getep_sub_tpc_recen[1][sub]     = (TProfile2D*)mRecenteringInputFile -> Get(Form("p_mq2y_tpc_AB_sub_%d",sub + 1));
      } // event plane Psi histograms
      getep_full_recen[0]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq1x_epd_ABCD_full");
      getep_full_recen[1]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq1y_epd_ABCD_full");
      
       // sub w/ eta weight 
      getep_full_wt_recen[0]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq1x_epd_ABCD_wt_full");
      getep_full_wt_recen[1]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq1y_epd_ABCD_wt_full");

       // tpc 
      getep_full_tpc_recen[0]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq2x_tpc_AB_full");
      getep_full_tpc_recen[1]     = (TProfile2D*)mRecenteringInputFile -> Get("p_mq2y_tpc_AB_full");
    cout << "Recenterpar successfully implemented" << endl;
    }
    //cout << (float)getep_full_recen[0]->GetBinContent(1, 1) << endl;
    //TFile *mShiftInputFile = TFile::Open("/star/data01/pwg/dchen/Ana/3GeV_FXT_2018/recenter/3gev_recenter.root");
    //TFile *mShiftInputFile = TFile::Open("/star/data01/pwg/dchen/Ana/19p6GeV/recenter/19gev_recenter.root"); // pwg
    TFile *mShiftInputFile = TFile::Open("/star/data01/pwg/dchen/Ana/19p6GeV/shift/19gev_shiftpar.root"); // pwg
    //TFile *mShiftInputFile = TFile::Open("/star/u/dchen/ana/19gev_2019/shift/test_shiftpar_shift.root"); // more test
    if (mShiftInputFile->IsZombie()) {
    	std::cout << "Error opening file with shift correction histograms" << std::endl;
    	std::cout << "No shifting in this analysis." << std::endl;
        for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
	{
          p_sub_ep_shiftpar_sin[sub] = 0;
          p_sub_ep_shiftpar_cos[sub] = 0;
          p_sub_ep_wt_shiftpar_sin[sub] = 0;
          p_sub_ep_wt_shiftpar_cos[sub] = 0;
          p_sub_ep_tpc_shiftpar_sin[sub] = 0;
          p_sub_ep_tpc_shiftpar_cos[sub] = 0;
	}
        p_full_ep_shiftpar_sin = 0;
        p_full_ep_shiftpar_cos = 0;
        p_full_ep_wt_shiftpar_sin = 0;
        p_full_ep_wt_shiftpar_cos = 0;
        p_full_ep_tpc_shiftpar_sin = 0;
        p_full_ep_tpc_shiftpar_cos = 0;
    }
    else{
      for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
      {
        p_sub_ep_shiftpar_sin[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_sub_%d_sin", sub + 1));
        p_sub_ep_shiftpar_cos[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_sub_%d_cos", sub + 1));
        p_sub_ep_wt_shiftpar_sin[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_wt_sub_%d_sin", sub + 1));
        p_sub_ep_wt_shiftpar_cos[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_wt_sub_%d_cos", sub + 1));
        p_sub_ep_tpc_shiftpar_sin[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_tpc_AB_sub_%d_sin", sub + 1));
        p_sub_ep_tpc_shiftpar_cos[sub]     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_tpc_AB_sub_%d_cos", sub + 1));
      }
      p_full_ep_shiftpar_sin     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_full_sin"));
      p_full_ep_shiftpar_cos     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_full_cos"));
      p_full_ep_wt_shiftpar_sin     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_wt_full_sin"));
      p_full_ep_wt_shiftpar_cos     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_epd_ABCD_wt_full_cos"));
      p_full_ep_tpc_shiftpar_sin     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_tpc_AB_full_sin"));
      p_full_ep_tpc_shiftpar_cos     = (TProfile3D*)mShiftInputFile -> Get(Form("EPshiftpar_tpc_AB_full_cos"));
      std::cout << "Shiftpar successfully implemented" << std::endl;
    }

    return 0;

}

//----------------------------------------------------------------------------- 
Int_t ParentParticle::Finish() 
{


    if(mOutPut_Phi != "")
    {
        mFile_Phi->cd();
        //mHistManager->WriteQAPID();
        //mEpProManager->WriteEP();
        mFlowV0->WritePhiMass2();
        mFile_Phi->Close();
    }
    return kStOK;
}

//----------------------------------------------------------------------------- 
void ParentParticle::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t ParentParticle::Make() 
{
    if(!mPicoDstMaker) 
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) 
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }
    //cout<<"test pp 1 "<<endl;
    const Int_t nTracks = mPicoDst->numberOfTracks();
    Int_t TrkMult = 0;
    for(Int_t i = 0; i < nTracks; i++) // TrkMult
    {
        StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
        if(!track->isPrimary()) continue; // Only Primary Tracks
        //mHistManager->FillTrackQA(track,(TVector3)mPicoEvent->primaryVertex());
        //StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
	TrkMult ++;
    } // TrkMult

    // RefMult
    Int_t runId = mPicoEvent->runId();

    Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    //Float_t vx = mPicoEvent->primaryVertex().X();
    //Float_t vy = mPicoEvent->primaryVertex().Y();

    //Float_t vzvpd = mPicoEvent->vzVpd();
    //Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    //Int_t nMatchedToF = mPicoEvent->nBTOFMatch();
    Float_t zdcX = mPicoEvent->ZDCx();
    
    //mHistManager->FillEventQA(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
    //mHistManager->FillEventCut(0);

    // vz sign
    Int_t vz_sign;
    if(vz > 0.0)
    {
        vz_sign = 0;
    }
    else
    {
        vz_sign = 1;
    }
    // runIndex
    const int runIndex = GetRunIndex(runId);
    // event plane IClasses
    //Start with clean events
    IEvent * theEvent = new IEvent;   
    IEvent * theEvent_wt = new IEvent;   
    IEvent * theEvent_tpc = new IEvent;   
    IEvent subEvents[_numSubEvents];
    IEvent subEvents_wt[_numSubEvents];
    IEvent subEvents_tpc[_numSubEvents];
    theEvent->ClearEvent();
    theEvent_wt->ClearEvent();
    theEvent_tpc->ClearEvent();
    for (int i = 0; i < _numSubEvents; i++)
    {
    	subEvents[i].ClearEvent();
    	subEvents_wt[i].ClearEvent();
    	subEvents_tpc[i].ClearEvent();
    }	

    // Event Cut
    if(mCutManager->passEventCut(mPicoDst)) // event cut
    {  
        //cout << " begin of a good event  " << endl;
	
        //mHistManager->FillEventQaCut(mPicoEvent->primaryVertex(),refMult,TOF_Mul,TrkMult);
        //mHistManager->FillEventCut(1);
        mRefMultCorr->init(runId);
        mRefMultCorr->initEvent(refMult, vz, zdcX);
        //const Int_t cent16 = mRefMultCorr->getCentralityBin16();
        const Int_t cent9 = mRefMultCorr->getCentralityBin9();
        //const Double_t reweight = mRefMultCorr->getWeight();
	//std::cout << "refMult: " << refMult << std::endl;
	//std::cout << "TrkMult: " << TrkMult << std::endl;
        //const int cent16 = mCutManager->getCentrality(TrkMult);
	
	//std::cout << "cent16: " << cent16 << std::endl;
	//std::cout << "cent9: " << cent9 << std::endl;
        //const double reweight = 1.0;
        //if(cent16 >  15 || cent16 < 0) return 0;
        if(cent9 > 8  || cent9 < 0) return 0;
        //mHistManager->FillEventCent(cent9);
        //mHistManager->FillEventCut(2);
	
	/// qapid
        //const Int_t nToFMatched = mCutManager->getMatchedToF();
        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        //float mField = mPicoEvent->bField();
        //Int_t N_pp = 0, N_pm = 0, N_kp = 0, N_km = 0, N_pr = 0;

	/// eventplane
	// TPC eventplane
        //double qx_TPC_A=0.0, qy_TPC_A=0.0;
        //double qx_TPC_B=0.0, qy_TPC_B=0.0;
	//
	// parentparticle
        float nMatchedToF=0, nN_prim=0, nN_non_prim=0;
        nMatchedToF = mPicoEvent->nBTOFMatch();

        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            if(!track->isPrimary()) continue; // Only Primary Tracks
            //mHistManager->FillTrackCut(0);
            if(!mCutManager->passTrackBasic(track)) continue;
            //mHistManager->FillTrackCut(1);
            //mHistManager->FillTrackPhysics(track );
    	    //std::cout << "pass track basic cut" << std::endl;
		
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            //Short_t  s_charge = track->charge();
            Float_t dca=track->gDCA(mVertexPos).Mag();
            if(dca > 3){// non primary track
                nN_non_prim++;
            }
            else{
                nN_prim++;
            }
            if(mCutManager->passTrackEP(track,dca)) // track cut   //shaowei
            {
                //cout<<"test tt 1 "<<endl;
                // calculate Q Vector after recentering for full event and eta sub event
                if(mFlowCorrection->passTrackFull(track))
                {
                    //cout<<"test tt 2 "<<endl;
                    mFlowCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
                    //cout<<"test tt 3 "<<endl;
                }
                for(Int_t j = 0; j < 4; j++) // eta_gap loop
                {
                    if(mFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
                    {
                    //cout<<"test tt 4 "<<endl;
                        mFlowCorrection->addTrack_East(track,cent9,runIndex,vz_sign,j);
                    //cout<<"add track east done! "<<endl;
                    }
                    if(mFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
                    {
                    //cout<<"test tt 5 "<<endl;
                        mFlowCorrection->addTrack_West(track,cent9,runIndex,vz_sign,j);
                    //cout<<"add track west done! "<<endl;
                    }
                }
            }
	    //if(mCutManager->isTofTrack(mPicoDst,track)) mHistManager->FillTrackTof(mPicoDst,track);
	    /*if(mCutManager->isProton(track))
	    {
	   	mHistManager->FillProton(mPicoDst,track,configs.y_mid); 
		N_pr++;
	    }
	    if(mCutManager->isKaon(mPicoDst,track))
	    {
	   	mHistManager->FillKaon(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_kp++;
		else N_km++;
	    }
	    if(mCutManager->isPion(mPicoDst,track))
	    {
	   	mHistManager->FillPion(mPicoDst,track,configs.y_mid); 
		if(s_charge >0) N_pp++;
		else N_pm++;
	    }
	    */
	    //TPC EP
            if(!mCutManager->passTrackEP(track,dca)) continue;
	    Double_t eta = track->pMom().Eta();
	    Double_t phi = track->pMom().Phi();
	    Double_t pt = track->pMom().Perp();
	    IEventPlane eventPlane(phi, pt);
	    eventPlane.SetEta(eta);
	    theEvent_tpc->AddEPParticle(eventPlane);
            /*if(eta > -2.0 && eta < -1.1)    // sub event A eta(-2,-1.25)   // mid-rapidity is 1.06  -2.0-1.4, -1.05-0.65
            {
	   	mEpProManager->FillTpcAQvec(cent16,runIndex, track); 
            }

            if(eta > -1.0 && eta < 0.)    // sub event B eta(-1.15,0) // -1.3-0.685,  -0.65-0
            {
	   	mEpProManager->FillTpcBQvec(cent16,runIndex,track); 
            }*/
        } // track loop

	/// eventplane
	// EPD information
        Int_t nepdHits = mPicoDst->numberOfEpdHits();
        StPicoEpdHit *epdHit;
        TVector3 StraightLine_center;
        TVector3 StraightLine_random;
        double phi_epd_center = {0.0};
        double phi_epd_random = {0.0};
        double eta_epd_center = {0.0};
        double eta_epd_random = {0.0};
        double mip;
        double TileWeight           = {0};

        //double Qx_EPD[7]={0};
        //double Qy_EPD[7]={0};
        //double weight_EPD[7]={0};
        for(Int_t iHit=0; iHit<nepdHits; iHit++){ // EPD loop
            epdHit = mPicoDst->epdHit(iHit);
            mip = epdHit->nMIP();
            int iring = epdHit->row() -1;//(1~16)-1 -> 0-15

            if( !epdHit) continue;
            if( !epdHit->isGood())
	    { 
            	std::cout << "note good epd hit "  << std::endl;
		    continue;
	    }
            int position = epdHit->position();
            int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner mose
            //std::cout << "ringgroup = " << ringgroup << std::endl;
            if(ringgroup == -1) continue;
            //std::cout << "epdHit id = " << epdHit->id() << std::endl;
            //if(epdHit->id() > 0 ) continue; // only East side for FXT 

            StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - mPicoEvent->primaryVertex();  // the collision is from midille to side of the TPC
            StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - mPicoEvent->primaryVertex();
            //StraightLine_center = mEpdGeom->TileCenter(epdHit->id())       ;
            //StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id());

            phi_epd_center = StraightLine_center.Phi();
            eta_epd_center = StraightLine_center.Eta();
            phi_epd_random = StraightLine_random.Phi();
            eta_epd_random = StraightLine_random.Eta();
	    mEpProManager->FillEpdQa(iring, phi_epd_center, eta_epd_center, phi_epd_random, eta_epd_random);
            if(mip < configs.epd_threshold) continue;
            TileWeight = (mip > configs.epd_max_weight) ? configs.epd_max_weight : mip;
            //std::cout << "tileweight  = " << TileWeight << std::endl;

	    mEpProManager->FillEpdMip(eta_epd_center,position,TileWeight);
	    Double_t eta_wt = mEpProManager->GetEtaWeight(cent9,eta_epd_center);
	    //std::cout<< "centrality : "<< cent9 << "eta: " << eta_epd_center << "eta weight: " << eta_wt << std::endl;

	    IEventPlane eventPlane(phi_epd_center, TileWeight);
	    IEventPlane eventPlane_wt(phi_epd_center, TileWeight*eta_wt);
	    //eventPlane.SetTileID(epdHit->id());
	    //std::cout << eta_epd_center << std::endl;
	    eventPlane.SetEta(eta_epd_center);
	    eventPlane_wt.SetEta(eta_epd_center);
	    theEvent->AddEPParticle(eventPlane);
	    theEvent_wt->AddEPParticle(eventPlane_wt);
            /*for(int i=0; i<4; i++){
                if(i == ringgroup){
                    Qx_EPD[i] += TileWeight * cos(1.0*phi_epd_center);
                    Qy_EPD[i] += TileWeight * sin(1.0*phi_epd_center);
                    weight_EPD[i] += TileWeight;
                }
            }
            if(ringgroup == 0 || ringgroup == 1)
            {
                Qx_EPD[4] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[4] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[4] += TileWeight;
            }
            if(ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[5] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[5] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[5] += TileWeight;
            }
            if(ringgroup == 0 || ringgroup == 1 || ringgroup == 2 || ringgroup == 3)
            {
                Qx_EPD[6] += TileWeight * cos(1.0*phi_epd_center);
                Qy_EPD[6] += TileWeight * sin(1.0*phi_epd_center);
                weight_EPD[6] += TileWeight;
            }*/
        } // EPD loop
    	/*for(int i=0; i<7; i++){
    	    if(weight_EPD[i] == 0. || Qx_EPD[i]== 0.0 || Qy_EPD[i] == 0.0) return 0;
    	    Qx_EPD[i] /= weight_EPD[i];
    	    Qy_EPD[i] /= weight_EPD[i];
    	}
    	for(int i=0; i<7; i++){
	    mEpProManager->FillEpdQvec(i,cent16,runIndex,Qx_EPD[i],Qy_EPD[i]);
	}*/
	//double psi_ievt = theEvent->GetEventPsi(1);
	//std::cout << "ievent psi = " << psi_ievt << std::endl;
	//mHistManager->FillPIDMult(N_pp , N_pm , N_kp , N_km , N_pr ); /// qapid
	for (int sub = 0; sub < _numSubEvents; sub++)
	{
	  //subEvents[sub] = theEvent->GetSubEvent(subEventModes[sub], subEventParams[sub][0] - COMrapidity, subEventParams[sub][1] - COMrapidity);
	  subEvents[sub] = theEvent->GetSubEvent(_subEventModes[sub], _subEventParams[sub][0], _subEventParams[sub][1]);
	  subEvents_wt[sub] = theEvent_wt->GetSubEvent(_subEventModes[sub], _subEventParams[sub][0], _subEventParams[sub][1]);
	  subEvents_tpc[sub] = theEvent_tpc->GetSubEvent(_subEventModes[sub], _subEventParams_tpc[sub][0], _subEventParams_tpc[sub][1]);
	}
	//std::cout << "Size of wt TPC-East sub EP = " << subEvents_tpc[0].GetEPParticles().size() << std::endl;
	//std::cout << "Size of wt TPC-West sub EP = " << subEvents_tpc[1].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  EPD-C sub EP = " << subEvents[0].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-A sub EP = " << subEvents[1].GetEPParticles().size() << std::endl;
	//std::cout << "Size of  TPC-B sub EP = " << subEvents[2].GetEPParticles().size() << std::endl << std::endl;
	if (subEvents[0].GetEPParticles().size() < 2
	//|| subEvents[1].GetEPParticles().size() < 2
	|| subEvents[1].GetEPParticles().size() < 2)
		{return 0;}	
	
	float subQx[_numSubEvents];
	float subQy[_numSubEvents];
	float subQx_wt[_numSubEvents];
	float subQy_wt[_numSubEvents];
	float subQx_tpc[_numSubEvents];
	float subQy_tpc[_numSubEvents];
	mEpProManager->FillPsiRawSubs(subEvents[0].GetEventPsi(1),subEvents[1].GetEventPsi(1));
	mEpProManager->FillPsiRawSubs_wt(subEvents_wt[0].GetEventPsi(1),subEvents_wt[1].GetEventPsi(1));
	mEpProManager->FillPsiRawSubs_tpc(subEvents_tpc[0].GetEventPsi(2),subEvents_tpc[1].GetEventPsi(2));
	
	mEpProManager->FillPsiRawFull(theEvent->GetEventPsi(1));
	mEpProManager->FillPsiRawFull_wt(theEvent_wt->GetEventPsi(1));
	mEpProManager->FillPsiRawFull_tpc(theEvent_tpc->GetEventPsi(2));
	
	
	mEpProManager->FillSubEpQvec_full(cent9, runIndex, theEvent->GetQx(1), theEvent->GetQy(1));
	mEpProManager->FillSubEpQvec_wt_full(cent9, runIndex, theEvent_wt->GetQx(1), theEvent_wt->GetQy(1));
	mEpProManager->FillSubEpQvec_tpc_full(cent9, runIndex, theEvent_tpc->GetQx(2), theEvent_tpc->GetQy(2));
	for (int sub = 0; sub < _numSubEvents; sub++) // raw event plane, store recenter parameter
	{
	  mEpProManager->FillPsiRaw(sub, subEvents[sub].GetEventPsi(1));
	  mEpProManager->FillPsiRaw_wt(sub, subEvents_wt[sub].GetEventPsi(1));
	  mEpProManager->FillPsiRaw_tpc(sub, subEvents_tpc[sub].GetEventPsi(2));
     	  subQx[sub] = subEvents[sub].GetQx(1);
     	  subQy[sub] = subEvents[sub].GetQy(1);
     	  subQx_wt[sub] = subEvents_wt[sub].GetQx(1);
     	  subQy_wt[sub] = subEvents_wt[sub].GetQy(1);
     	  subQx_tpc[sub] = subEvents_tpc[sub].GetQx(2);
     	  subQy_tpc[sub] = subEvents_tpc[sub].GetQy(2);
	  mEpProManager->FillSubEpQvec(sub, cent9, runIndex, subQx[sub], subQy[sub]);
	  mEpProManager->FillSubEpQvec_wt(sub, cent9, runIndex, subQx_wt[sub], subQy_wt[sub]);
	  mEpProManager->FillSubEpQvec_tpc(sub, cent9, runIndex, subQx_tpc[sub], subQy_tpc[sub]);
	} // raw event plane, store recenter parameter
	
	// EPD ABCD
	for (int i = 0; i < _numSubEvents; i++)
	{
	  if (getep_sub_recen[0][i] != 0 && getep_sub_recen[1][i] != 0) // sub recenter event plane, store shift parameter
	     {
	       float subEventQx = (float)getep_sub_recen[0][i]->GetBinContent(cent9+1, runIndex+1);
	       float subEventQy = (float)getep_sub_recen[1][i]->GetBinContent(cent9+1, runIndex+1);
	       subEvents[i].SetQCenter(subEventQx, subEventQy);
    	        //std::cout << subEventQx  << " subEventQx " << std::endl;
    	        //std::cout << subEventQy  << " subEventQy " << std::endl;
	     }
	} // sub recenter event plane, store shift parameter
	if (getep_full_recen[0] != 0 && getep_full_recen[1] != 0) // full recenter event plane, store shift parameter
        {
	       float fullEventQx = (float)getep_full_recen[0]->GetBinContent(cent9+1, runIndex+1);
	       float fullEventQy = (float)getep_full_recen[1]->GetBinContent(cent9+1, runIndex+1);
	       theEvent->SetQCenter(fullEventQx, fullEventQy);
        } // full recenter event plane, store shift parameter


	float subQxRec[_numSubEvents];
	float subQyRec[_numSubEvents];
	float psiShift[_numSubEvents];
	for (int i = 0; i < _numSubEvents; i++)
	{ // recenter event plane, store shift parameter
        	if(subEvents[i].GetEventPsi(1) == -99){continue;} //Empty subevent
	        mEpProManager->FillPsiRec(i, subEvents[i].GetEventPsi(1));
     	  	subQxRec[i] = subEvents[i].GetQx(1);
     	  	subQyRec[i] = subEvents[i].GetQy(1);
		// implement shift correction
		float f_Psi_shifted = subEvents[i].GetEventPsi(1);
	        mEpProManager->FillSubEpQvecRec(i, cent9, runIndex, subQxRec[i], subQyRec[i]);
                for(int iorder=0; iorder<20; iorder++)
                {
                    // first order shift parameter
	            mEpProManager->FillSubEpShiftpar(i, iorder, cent9, runIndex, subEvents[i].GetEventPsi(1));
                    // input first order shift parameter
		    // implement shift correction
	       	    float tmp = (float)(configs.order_m*(iorder+1));
	       	    float sinAve = p_sub_ep_shiftpar_sin[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       	    float cosAve = p_sub_ep_shiftpar_cos[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       	    f_Psi_shifted +=
	       	    2.0*(cosAve*TMath::Sin(tmp* subEvents[i].GetEventPsi(1)) - sinAve*TMath::Cos(tmp* subEvents[i].GetEventPsi(1)))/tmp;
                }
		double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_m;
		double HalfWrapAround = AngleWrapAround/2;
		while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
		while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	        mEpProManager->FillPsiShift(i, f_Psi_shifted);
		psiShift[i] = f_Psi_shifted;
	} // recenter event plane, store shift parameter
	mEpProManager->FillPsiRecSubs(subEvents[0].GetEventPsi(1),subEvents[1].GetEventPsi(1));
	mEpProManager->FillPsiShiftSubs(psiShift[0], psiShift[1]);
	// Fill resolution plots
	mEpProManager->FillPsi1ResolutionEpd(cent9, psiShift[0], psiShift[1]);
	
	// full event plane	
        if(theEvent->GetEventPsi(1) != -99)
	{ // recenter the full event plane, shift parameter
	   mEpProManager->FillPsiRecFull(theEvent->GetEventPsi(1));
	   float f_Psi_shifted = theEvent->GetEventPsi(1);
           for(int iorder=0; iorder<20; iorder++)
           {
               // first order shift parameter
	       mEpProManager->FillFullEpShiftpar(iorder, cent9, runIndex, theEvent->GetEventPsi(1));
               // input first order shift parameter
	       // implement shift correction
	       float tmp = (float)(configs.order_m*(iorder+1));
	       float sinAve = p_full_ep_shiftpar_sin->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       float cosAve = p_full_ep_shiftpar_cos->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       f_Psi_shifted +=
	       2.0*(cosAve*TMath::Sin(tmp* theEvent->GetEventPsi(1)) - sinAve*TMath::Cos(tmp* theEvent->GetEventPsi(1)))/tmp;
           }
	   double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_m;
	   double HalfWrapAround = AngleWrapAround/2;
	   while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
	   while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	   mEpProManager->FillPsiShiftFull( f_Psi_shifted);
        } // recenter the full event plane, shift parameter

	// EPD ABCD w/ wt
	for (int i = 0; i < _numSubEvents; i++)
	{
	  if (getep_sub_wt_recen[0][i] != 0 && getep_sub_wt_recen[1][i] != 0) // sub_wt recenter event plane, store shift parameter
	     {
	       float subEventQx = (float)getep_sub_wt_recen[0][i]->GetBinContent(cent9+1, runIndex+1);
	       float subEventQy = (float)getep_sub_wt_recen[1][i]->GetBinContent(cent9+1, runIndex+1);
	       subEvents_wt[i].SetQCenter(subEventQx, subEventQy);
	     }
	} // sub recenter event plane, store shift parameter
	if (getep_full_wt_recen[0] != 0 && getep_full_wt_recen[1] != 0) // full recenter event plane, store shift parameter
        {
	       float fullEventQx = (float)getep_full_wt_recen[0]->GetBinContent(cent9+1, runIndex+1);
	       float fullEventQy = (float)getep_full_wt_recen[1]->GetBinContent(cent9+1, runIndex+1);
	       theEvent_wt->SetQCenter(fullEventQx, fullEventQy);
        } // full recenter event plane, store shift parameter

	//float subQxRec_wt[_numSubEvents];
	//float subQyRec_wt[_numSubEvents];
	float psiShift_wt[_numSubEvents];
	for (int i = 0; i < _numSubEvents; i++)
	{ // recenter event plane, store shift parameter

        	if(subEvents_wt[i].GetEventPsi(1) == -99){continue;} //Empty subevent
	        mEpProManager->FillPsiRec_wt(i, subEvents_wt[i].GetEventPsi(1));
     	  	//subQxRec_wt[i] = subEvents_wt[i].GetQx(1);
     	  	//subQyRec_wt[i] = subEvents_wt[i].GetQy(1);
		// implement shift correction
		float f_Psi_shifted = subEvents_wt[i].GetEventPsi(1);
                for(int iorder=0; iorder<20; iorder++)
                {
                    // first order shift parameter
	            mEpProManager->FillSubEpShiftpar_wt(i, iorder, cent9, runIndex, subEvents_wt[i].GetEventPsi(1));
                    // input first order shift parameter
		    // implement shift correction
	       	    float tmp = (float)(configs.order_m*(iorder+1));
	       	    float sinAve = p_sub_ep_wt_shiftpar_sin[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       	    float cosAve = p_sub_ep_wt_shiftpar_cos[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       	    f_Psi_shifted +=
	       	    2.0*(cosAve*TMath::Sin(tmp* subEvents_wt[i].GetEventPsi(1)) - sinAve*TMath::Cos(tmp* subEvents_wt[i].GetEventPsi(1)))/tmp;
                }
		double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_m;
		double HalfWrapAround = AngleWrapAround/2;
		while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
		while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	        mEpProManager->FillPsiShift_wt(i, f_Psi_shifted);
		psiShift_wt[i] = f_Psi_shifted;
	} // recenter event plane, store shift parameter
	mEpProManager->FillPsiRecSubs_wt(subEvents_wt[0].GetEventPsi(1),subEvents_wt[1].GetEventPsi(1));
	mEpProManager->FillPsiShiftSubs_wt(psiShift_wt[0], psiShift_wt[1]);
	// Fill resolution plots
	mEpProManager->FillPsi1ResolutionEpdWt(cent9, psiShift_wt[0], psiShift_wt[1]);
	
	// full event plane	
        if(theEvent_wt->GetEventPsi(1) != -99)
	{ // recenter the full event plane, shift parameter
	   mEpProManager->FillPsiRecFull_wt(theEvent_wt->GetEventPsi(1));
	   float f_Psi_shifted = theEvent_wt->GetEventPsi(1);
           for(int iorder=0; iorder<20; iorder++)
           {
               // first order shift parameter
	       mEpProManager->FillFullEpShiftpar_wt(iorder, cent9, runIndex, theEvent_wt->GetEventPsi(1));
               // input first order shift parameter
	       // implement shift correction
	       float tmp = (float)(configs.order_m*(iorder+1));
	       float sinAve = p_full_ep_wt_shiftpar_sin->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       float cosAve = p_full_ep_wt_shiftpar_cos->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       f_Psi_shifted +=
	       2.0*(cosAve*TMath::Sin(tmp* theEvent->GetEventPsi(1)) - sinAve*TMath::Cos(tmp* theEvent->GetEventPsi(1)))/tmp;
	   }
	   double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_m;
	   double HalfWrapAround = AngleWrapAround/2;
	   while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
	   while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	   mEpProManager->FillPsiShiftFull_wt( f_Psi_shifted);
        } // recenter the full event plane, shift parameter

        // TPC AB recenter
	for (int i = 0; i < _numSubEvents; i++)
	{
	  if (getep_sub_tpc_recen[0][i] != 0 && getep_sub_tpc_recen[1][i] != 0) // sub_tpc recenter event plane, store shift parameter
	     {
	       float subEventQx = (float)getep_sub_tpc_recen[0][i]->GetBinContent(cent9+1, runIndex+1);
	       float subEventQy = (float)getep_sub_tpc_recen[1][i]->GetBinContent(cent9+1, runIndex+1);
		//cout << i <<  " " << subEventQx << " Qx "<< subEventQy << " Qy " << "cent "<< cent9  << "runIndex: "<< runIndex << "bin content: "<<  endl;
	       subEvents_tpc[i].SetQCenter(subEventQx, subEventQy);
	     }
	} // sub recenter event plane, store shift parameter
	if (getep_full_tpc_recen[0] != 0 && getep_full_tpc_recen[1] != 0) // full recenter event plane, store shift parameter
        {
	       float fullEventQx = (float)getep_full_tpc_recen[0]->GetBinContent(cent9+1, runIndex+1);
	       float fullEventQy = (float)getep_full_tpc_recen[1]->GetBinContent(cent9+1, runIndex+1);
	       theEvent_tpc->SetQCenter(fullEventQx, fullEventQy);
        } // full recenter event plane, store shift parameter

	//float subQxRec_tpc[_numSubEvents];
	//float subQyRec_tpc[_numSubEvents];
	float psiShift_tpc[_numSubEvents];
	for (int i = 0; i < _numSubEvents; i++)
	{ // recenter event plane, store shift parameter
        	if(subEvents_tpc[i].GetEventPsi(2) == -99){continue;} //Empty subevent
		//cout << i << subEvents_tpc[i].GetEventPsi(2) << " psi2 tpc sub" << endl;
	        mEpProManager->FillPsiRec_tpc(i, subEvents_tpc[i].GetEventPsi(2));
     	  	//subQxRec_tpc[i] = subEvents_tpc[i].GetQx(2);
     	  	//subQyRec_tpc[i] = subEvents_tpc[i].GetQy(2);
		// implement shift correction
		float f_Psi_shifted = subEvents_tpc[i].GetEventPsi(2);
                for(int iorder=0; iorder<20; iorder++)
                {
                    // first order shift parameter
	            mEpProManager->FillSubEpShiftpar_tpc(i, iorder, cent9, runIndex, subEvents_tpc[i].GetEventPsi(2));
                    // input first order shift parameter
		    // implement shift correction
	       	    float tmp = (float)(configs.order_n*(iorder+1));
	       	    float sinAve = p_sub_ep_tpc_shiftpar_sin[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       	    float cosAve = p_sub_ep_tpc_shiftpar_cos[i]->GetBinContent(cent9+1, iorder+1, runIndex+1);
		//cout << i <<  " " << p_sub_ep_tpc_shiftpar_cos[i] << " p_cos[i] " << "cent "<< cent9 << "iorder: "<< iorder << "runIndex: "<< runIndex << "bin content: "<<p_sub_ep_tpc_shiftpar_cos[i]->GetBinContent(cent9+1, iorder+1, runIndex+1)<<  endl;
	       	    f_Psi_shifted +=
	       	    2.0*(cosAve*TMath::Sin(tmp* subEvents_tpc[i].GetEventPsi(2)) - sinAve*TMath::Cos(tmp* subEvents_tpc[i].GetEventPsi(2)))/tmp;
                }
		double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_n;
		double HalfWrapAround = AngleWrapAround/2;
		while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
		while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	        mEpProManager->FillPsiShift_tpc(i, f_Psi_shifted);
		psiShift_tpc[i] = f_Psi_shifted;
		//cout << i <<  " " << f_Psi_shifted << " psi2 tpc sub" << endl;
	} // recenter event plane, store shift parameter
	mEpProManager->FillPsiRecSubs_tpc(subEvents_tpc[0].GetEventPsi(2),subEvents_tpc[1].GetEventPsi(2));
	mEpProManager->FillPsiShiftSubs_tpc(psiShift_tpc[0], psiShift_tpc[1]);
	// Fill resolution plots
	mEpProManager->FillPsi2ResolutionTpc(cent9, psiShift_tpc[0], psiShift_tpc[1]);
	//cout << "tpc east = " << psiShift_tpc[0] << " tpc west = "<<psiShift_tpc[1] << endl;

	// full event plane	
        if(theEvent_tpc->GetEventPsi(2) != -99)
	{ // recenter the full event plane, shift parameter
	   mEpProManager->FillPsiRecFull_tpc(theEvent_tpc->GetEventPsi(2));
	   float f_Psi_shifted = theEvent_tpc->GetEventPsi(2);
           for(int iorder=0; iorder<20; iorder++)
           {
               // second order shift parameter
	       mEpProManager->FillFullEpShiftpar_tpc(iorder, cent9, runIndex, theEvent_tpc->GetEventPsi(2));
	       // implement shift correction
	       float tmp = (float)(configs.order_n*(iorder+1));
	       float sinAve = p_full_ep_tpc_shiftpar_sin->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       float cosAve = p_full_ep_tpc_shiftpar_cos->GetBinContent(cent9+1, iorder+1, runIndex+1);
	       f_Psi_shifted +=
	       2.0*(cosAve*TMath::Sin(tmp* theEvent->GetEventPsi(2)) - sinAve*TMath::Cos(tmp* theEvent->GetEventPsi(2)))/tmp;
	   }
	   double AngleWrapAround = 2.0*TMath::Pi()/(double)configs.order_n;
	   double HalfWrapAround = AngleWrapAround/2;
	   while (f_Psi_shifted < -1*HalfWrapAround) f_Psi_shifted += AngleWrapAround;
	   while (f_Psi_shifted > HalfWrapAround) f_Psi_shifted -= AngleWrapAround;
	   mEpProManager->FillPsiShiftFull_tpc( f_Psi_shifted);
        } // recenter the full event plane, shift parameter
        Float_t Psi2_East;
        TVector2 Q1Full[ConstManager::EtaGap_total];
        TVector2 Q1EpdWt;
        Q1EpdWt.Set(theEvent_wt->GetQx(1),theEvent_wt->GetQy(1)); // Epd wt full eventplane
        TVector2 Q2East[ConstManager::EtaGap_total], Q2West[ConstManager::EtaGap_total];
        TVector2 Q3East[ConstManager::EtaGap_total], Q3West[ConstManager::EtaGap_total];
        Int_t NumTrackEast[ConstManager::EtaGap_total], NumTrackWest[ConstManager::EtaGap_total]; 
        for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
        {
            Q1Full[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
            Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
            Q2West[j].Set(-999.9,-999.9);
            Q3East[j].Set(-999.9,-999.9);
            Q3West[j].Set(-999.9,-999.9);
            NumTrackEast[j] = 0;
            NumTrackWest[j] = 0;
            Q1Full[j]=Q1EpdWt; // 
	    //cout << "epd psi 1 wt qx: "<< theEvent_wt->GetQx(1) << "epd psi 1 wt qy: " << theEvent_wt->GetQy(1) << endl;
            if(mFlowCorrection->passTrackEtaNumCut(j))
            {
                // get QVector of sub event
                Q2East[j] = mFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q2West[j] = mFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q3East[j] = mFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q3West[j] = mFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
                NumTrackEast[j] = mFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
                NumTrackWest[j] = mFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
            }
        }
        //cout<<"event plane success! "<<endl;
     	/*for(int i=0; i<4; i++)	
	{
           Q1East[i] = subEvents_wt[0].GetQx(1);
           Q1West[i] = subEvents_wt[1].GetQy(1);
           Q2East[i] = subEvents_tpc[0].GetQx(2);
           Q2East[i] = subEvents_tpc[1].GetQy(2);
	}*/
        if(mFlowCorrection->passTrackEtaNumCut(0))
        {
            // Event Plane method
	    //cout << "begin track" << endl;
            Psi2_East = subEvents_tpc[0].GetEventPsi(2);//mFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);
            //Psi1_East = psiShift_wt[0];// mFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);

            // get N_prim, N_non_prim, N_Tof_match
            Int_t N_prim = nN_prim;
            Int_t N_non_prim = nN_non_prim;
            Int_t N_Tof_match = nMatchedToF;
            // pass the event information to StFlowV0
            mFlowV0->clearEvent();
            mFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

            // 1st sub event plane
            mFlowV0->passEventPlane1Full(Q1Full[0],Q1Full[1],Q1Full[2],Q1Full[3]);
            //mFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

            // 2nd sub event plane
            mFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
            mFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

            // 3rd sub event plane
            mFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
            mFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);

            // Number of Track in East and West part of TPC
            mFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
            mFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);

            //FlowV0->MixEvent_Phi(configs::Flag_ME,mPicoDst,cent9,vz,Psi2_East,reweight);
            //cout<<"begin phi meson reconstruction"<<endl;
	    //cout << configs::Flag_ME << " "<< mPicoDst << " "<< vz << " "<< " "<< Psi2_East << endl;
            mFlowV0->MixEvent_Phi(configs.Flag_ME,mPicoDst,cent9,vz,Psi2_East,1.0);
            //cout<<"finish phi meson reconstruction"<<endl;
        }
        mFlowCorrection->clear();

    //cout << " end of the good event " << endl;

    } // event cut

    return kStOK;
}

int ParentParticle::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<nrun; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}

/*int ParentParticle::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={4, 9,17,30,50,78, 116,170,205};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}*/
