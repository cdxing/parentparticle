#include "FlowV0.h"
//#include "StTriFlowConstants.h"
#include "StRoot/ConstManager/ConstManager.h" // dchen
//#include "StTriFlowCut.h"
#include "StRoot/CutManager/CutManager.h" // dchen
#include "StRoot/StPicoEvent/StPicoDst.h"   //shaowei
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"  // shaowei
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
//#include "StRoot/StV0Event/StV0Event.h"
//#include "StRoot/StV0TofCorrection/StV0TofCorrection.h"
#include <vector>
#include "TLorentzVector.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TTree.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "StLorentzVectorD.hh"

ClassImp(FlowV0)

    //------------------------------------------------------------------------------------------------------------------
//FlowV0::FlowV0(Int_t energy)
FlowV0::FlowV0(ConfigReader configs)
{
    mConfigs = configs;
}

FlowV0::~FlowV0()
{
}

//------------------------------------------------------------------------------------------------------------------

void FlowV0::InitPhi()
{
    // PID plots
    h_pT_kp = new TH1D("h_pT_kp","K^{+} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
    h_pT_km = new TH1D("h_pT_km","K^{-} p_{T}; p_{T} (GeV/c);Tracks",1000,0.0,5.0);
    h_eta_kp = new TH1D("h_eta_kp","K^{+} #eta;#eta;Tracks",500,-5.0,5.0);
    h_eta_km = new TH1D("h_eta_km","K^{-} #eta;#eta;Tracks",500,-5.0,5.0);
    h_dndy_kp = new TH1D("h_dndy_kp", "K^{+} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
    h_dndy_km = new TH1D("h_dndy_km", "K^{-} Raw Rapidity Spectrum;y;dN/dy",   80, -2, 2);
    h_phi_kp = new TH1D("h_phi_kp","K^{+} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
    h_phi_km = new TH1D("h_phi_km","K^{-} #phi (Radian);#phi;Tracks",1000,-1.5*TMath::Pi(),1.5*TMath::Pi());
    h_mult_kp = new TH1D("h_mult_kp","K^{#plus} track multiplicity;K^{+} Mult;Events",1001,-0.5,1000.5);
    h_mult_km = new TH1D("h_mult_km","K^{-} track multiplicity;K^{-} Mult;Events",1001,-0.5,1000.5);
    h2_dEdx_vs_qp_kp = new TH2D("h2_dEdx_vs_qp_kp", "K^{+} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
    h2_dEdx_vs_qp_km = new TH2D("h2_dEdx_vs_qp_km", "K^{-} dE/dx vs q|p|;q|p| (GeV);dE/dx (keV/cm)", 400, -2, 2, 500, 0, 10);
    h2_beta_vs_qp_kp = new TH2D("h2_beta_vs_qp_kp","K^{+} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
    h2_beta_vs_qp_km = new TH2D("h2_beta_vs_qp_km","K^{-} 1/#beta vs Momentum;q*|p| (GeV);1/#beta", 300, -3, 3, 300, 0.5, 3.5);
    h2_m2_vs_qp_kp = new TH2D("h2_m2_vs_qp_kp", "K^{+} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
    h2_m2_vs_qp_km = new TH2D("h2_m2_vs_qp_km", "K^{-} m^2 vs q*|p|;q*|p| (GeV);m^2 (GeV^2)",   400, -4, 4, 400, -0.1, 1.5);
    h2_pT_vs_yCM_kp = new TH2D("h2_pT_vs_yCM_kp", "K^{+};y-y_{mid};p_{T} (GeV/c)",   300, -2., 2., 300, 0, 3.0);
    h2_pT_vs_yCM_km = new TH2D("h2_pT_vs_yCM_km", "K^{-};y-y_{mid};p_{T} (GeV/c)",   300, -1.2, 1.2, 300, 0, 3.0);
    //mCutManager = new StTriFlowCut(mEnergy);
    mCutManager = new CutManager(mConfigs);
    TString HistName = "Mass2_pt";
    h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.98,1.08);

    for(Int_t cent = 0; cent < ConstManager::Bin_Centrality; cent++)
    {
        for(Int_t vz = 0; vz < ConstManager::Bin_VertexZ; vz++)
        {
            for(Int_t phi_psi = 0; phi_psi < ConstManager::Bin_Phi_Psi; phi_psi++)
            {
                mEventCounter2[cent][vz][phi_psi] = 0;
                clear_phi(cent,vz,phi_psi);
            }
        }
    }

    mXuPhiMesonEvent = new StAlexPhiMesonEvent();
    mTree_Phi = new TTree("XuPhiMesonEvent","XuPhiMesonEvent");
    mTree_Phi->Branch("phi_flow_branch","StAlexPhiMesonEvent",&mXuPhiMesonEvent);
    mTree_Phi->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------

void FlowV0::WritePhiMass2()
{
    h_mult_kp->Write();
    h_mult_km->Write();
    h_eta_kp->Write();
    h_phi_kp->Write();
    h_pT_kp->Write();
    h_dndy_kp->Write();
    h2_pT_vs_yCM_kp->Write();
    h2_dEdx_vs_qp_kp->Write();
    h2_beta_vs_qp_kp->Write();
    h2_m2_vs_qp_kp->Write();

    h_eta_km->Write();
    h_phi_km->Write();
    h_pT_km->Write();
    h_dndy_km->Write();
    h2_pT_vs_yCM_km->Write();
    h2_dEdx_vs_qp_km->Write();
    h2_beta_vs_qp_km->Write();
    h2_m2_vs_qp_km->Write();
    h_Mass2->Write();
    mTree_Phi->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void FlowV0::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
    mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
    //mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
    mReweight[cent9][Bin_vz][Bin_Psi2].clear();
    mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
    mRunId[cent9][Bin_vz][Bin_Psi2].clear();
    //mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
    //mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
    //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
    //mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
    //mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
    //mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();

    for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
    {
        mQ1Full[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
    }

    for(Int_t Bin_Event = 0; Bin_Event < ConstManager::Buffer_depth; Bin_Event++)
    {
        for(Int_t charge = 0; charge < 2; charge++)
        {
            MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
            mHelix_Kaon[key].clear();
            mMomentum[key].clear();
            mMass2[key].clear();
            mDca[key].clear();
            mNHitsFit[key].clear();
            mNSigmaKaon[key].clear();
        }
    }
    mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void FlowV0::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
    LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", Psi2 = " << Bin_Psi2 << endm;
    LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;

    LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    //LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    LOG_INFO << "Size of reweight       = " << mReweight[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    LOG_INFO << "---------------------------------------------------------------------------" << endm;

    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
        MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
        LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
        LOG_INFO << "Positive Particle:" << endm;
        LOG_INFO << "  Size of Helix_Kplus  = " << mHelix_Kaon[key].size() << endm;;
        LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
        LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
        LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
        LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
        LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;

        key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
        LOG_INFO << "Negative Particle:" << endm;
        LOG_INFO << "  Size of Helix_Kminus = " << mHelix_Kaon[key].size() << endm;
        LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
        LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
        LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
        LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
        LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;
        LOG_INFO << "---------------------------------------------------------------------------" << endm;
    }
}

//------------------------------------------------------------------------------------------------------------------

void FlowV0::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
    if(Flag_ME == 0) // same event
    {
            //cout << "Bin_Event = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endl;
        for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
        {
            //cout << "Bin_Event start  = " <<  endl;
            // event header
            //cout << "Bin_Event test 1 = " <<  endl;
	    //cout << "clear event success!" << endl;
            //cout << "Num tracks: "<< mXuPhiMesonEvent->getNumTracks() << endl;
	    //cout << "get numbers  success!" << endl;
            mXuPhiMesonEvent->clearTrackList();
            mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //cout << "Bin_Event test 2 = " <<  endl;
            mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //cout << "Bin_Event test 3 = " <<  endl;
            //mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setReweight(mReweight[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

            //cout << "Bin_Event test 0  " <<  endl;
            for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
            {
                // QVector
                mXuPhiMesonEvent->setQ1Full(mQ1Full[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                /*mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);*/
                // Number of Tracks
                mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
            }
            //cout << "Bin_Event test 1  " <<  endl;

            //mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            //mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

            // start to select phi candidate in a event
            MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
            MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

            TLorentzVector ltrackA, ltrackB;
            //cout <<"K+ number = " <<  mHelix_Kaon[key_plus].size() << endl;
            for(unsigned int n_kplus = 0; n_kplus < mHelix_Kaon[key_plus].size(); n_kplus++) // first track loop over K+ candidates
            {
                TVector3 p_vecA = mHelix_Kaon[key_plus][n_kplus].cat(mHelix_Kaon[key_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
                p_vecA *= mMomentum[key_plus][n_kplus];
                ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),ConstManager::mMassKaon);

                for(unsigned int n_kminus = 0; n_kminus < mHelix_Kaon[key_minus].size(); n_kminus++) // second track loop over K- candidates
                {
                    TVector3 p_vecB = mHelix_Kaon[key_minus][n_kminus].cat(mHelix_Kaon[key_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
                    p_vecB *= mMomentum[key_minus][n_kminus];
                    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),ConstManager::mMassKaon);

                    TLorentzVector trackAB      = ltrackA+ltrackB;
                    Double_t InvMassAB          = trackAB.M();
                    Double_t pt = trackAB.Perp();

                    // fill phi candidate into mTree_Phi
                    if(InvMassAB > ConstManager::mMassKaon*2 && InvMassAB < 1.05) 
                    {
                        mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                        mXuPhiMesonTrack->setMass2A(mMass2[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setMass2B(mMass2[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setDcaA(mDca[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setDcaB(mDca[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setTrackA(ltrackA); // K+
                        mXuPhiMesonTrack->setTrackB(ltrackB); // K-
                        mXuPhiMesonTrack->setFlagA(Bin_Event); // K+
                        mXuPhiMesonTrack->setFlagB(Bin_Event); // K-
                    }

                    // Fill histogram with InvMassAB information
                    h_Mass2->Fill(pt,InvMassAB);
                }
            }
        }
        mTree_Phi->Fill();
    }

    if(Flag_ME == 1) // mixed event
    {
        for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
        {
            MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
            MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
            for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
            {
                MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
                MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

                if(Bin_Event_A == 0 && Bin_Event_B == 1)
                {
                    Int_t Bin_Event = Bin_Event_A;
                    // event header
                    mXuPhiMesonEvent->clearTrackList();
                    mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setReweight(mReweight[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

                    for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
                    {
                        // QVector
                        mXuPhiMesonEvent->setQ1Full(mQ1Full[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        /*mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);*/

                        // Number of Tracks
                        mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                    }

                    //mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    //mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                }

                TLorentzVector ltrackA, ltrackB;

                // start to mix events
                // mix K+ candidates from A event with K- candidates from B event
                for(unsigned int n_kplus = 0; n_kplus < mHelix_Kaon[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event A
                {
                    TVector3 p_vecA(mHelix_Kaon[key_A_plus][n_kplus].cat(mHelix_Kaon[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
                    p_vecA *= mMomentum[key_A_plus][n_kplus];
                    ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),ConstManager::mMassKaon); // K+

                    for(unsigned int n_kminus = 0; n_kminus < mHelix_Kaon[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from event B
                    {
                        TVector3 p_vecB(mHelix_Kaon[key_B_minus][n_kminus].cat(mHelix_Kaon[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
                        p_vecB *= mMomentum[key_B_minus][n_kminus];
                        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),ConstManager::mMassKaon); // K-

                        TLorentzVector trackAB      = ltrackA+ltrackB;
                        Double_t InvMassAB          = trackAB.M();
                        Double_t pt = trackAB.Perp();

                        // fill phi candidate background into mTree_Phi
                        if(InvMassAB > ConstManager::mMassKaon*2 && InvMassAB < 1.05) 
                        {
                            mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                            mXuPhiMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setTrackA(ltrackA); // K+
                            mXuPhiMesonTrack->setTrackB(ltrackB); // K-
                            mXuPhiMesonTrack->setFlagA(Bin_Event_A); // K+
                            mXuPhiMesonTrack->setFlagB(Bin_Event_B); // K-
                        }

                        // Fill histogram with InvMassAB information
                        h_Mass2->Fill(pt,InvMassAB);
                    }
                }

                // mix K- candidates from A event with K+ candidates from B event
                for(unsigned int n_kminus = 0; n_kminus < mHelix_Kaon[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event A
                {
                    TVector3 p_vecA(mHelix_Kaon[key_A_minus][n_kminus].cat(mHelix_Kaon[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
                    p_vecA *= mMomentum[key_A_minus][n_kminus];
                    ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),ConstManager::mMassKaon); // K-

                    for(unsigned int n_kplus = 0; n_kplus < mHelix_Kaon[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from event B
                    {
                        TVector3 p_vecB(mHelix_Kaon[key_B_plus][n_kplus].cat(mHelix_Kaon[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
                        p_vecB *= mMomentum[key_B_plus][n_kplus];
                        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),ConstManager::mMassKaon); // K+

                        TLorentzVector trackAB      = ltrackA+ltrackB;
                        Double_t InvMassAB          = trackAB.M();
                        Double_t pt = trackAB.Perp();

                        // fill phi candidate background into mTree_Phi
                        if(InvMassAB > ConstManager::mMassKaon*2 && InvMassAB < 1.05) 
                        {
                            mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                            mXuPhiMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setTrackA(ltrackB); // K+
                            mXuPhiMesonTrack->setTrackB(ltrackA); // K-
                            mXuPhiMesonTrack->setFlagA(Bin_Event_B); // K+
                            mXuPhiMesonTrack->setFlagB(Bin_Event_A); // K-
                        }

                        // Fill histogram with InvMassAB information
                        h_Mass2->Fill(pt,InvMassAB);
                    }
                }
            }
        }
        mTree_Phi->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------


void FlowV0::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2, Float_t reweight)
{
	//cout << "begin MixEvent_Phi " << endl;
    StPicoEvent *event = (StPicoEvent*)pico->event();
    Int_t Bin_vz, Bin_Psi2;

    //Float_t vz_start = ConstManager::mVzMaxMap[event->energy()];
    Float_t vz_start = mConfigs.z_vtx_high; //   shaowei   same with 39 62GeV
    Float_t vz_bin = 2*vz_start/ConstManager::Bin_VertexZ;

    Float_t psi2_start = TMath::Pi()/2.0;
    Float_t psi2_bin = 2*psi2_start/ConstManager::Bin_Phi_Psi;
	//cout << "mix event v0 test 1" << endl;

    for(Int_t i = 0; i < ConstManager::Bin_VertexZ; i++)
    {
        if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
        {
            Bin_vz = i;
        }
    }
    for(Int_t i = 0; i < ConstManager::Bin_Phi_Psi; i++)
    {
        if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
        {
            Bin_Psi2 = i;
        }
    }
	//cout << "mix event v0 test 2" << endl;
	//cout << cent9 << " "<< Bin_vz << " " << Bin_Psi2 << endl;

    //mEventCounter2[cent9][Bin_vz][Bin_Psi2] =0;
    //cout <<"Bin_Event = "<< mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endl;
    Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];
	//cout << "mix event v0 test 3" << endl;

    const Double_t MAGFIELDFACTOR = kilogauss;
    const Int_t nTracks = pico->numberOfTracks();

    // store Enent Information
    mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector3>(event->primaryVertex()));
    //mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
    mReweight[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(reweight));
    mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
    mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
    //mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
    //mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
    //mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
    //mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
    //mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
    //mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
    mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
    //cout << "mix event v0 test 4" << endl;
    for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
    {
        mQ1Full[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector1Full[j]));
        mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
        mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
        mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
        mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
        mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
        mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
    }
	//cout << "mix event v0 test 5" << endl;

    // store Track Information
    TVector3 mVertexPos = event->primaryVertex();
    float mField = event->bField();
    //cout << "total track = " << nTracks << endl;
    //cout << "total tof track =  " << pico->numberOfBTofPidTraits() << endl;
    for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
    {
        StPicoTrack *track = pico->track(i);
        if(!track)
        {
            cout << "No PicoTrack Skip! " << endl;
            return;
        }
        StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
        Float_t dca=track->gDCA(mVertexPos).Mag();

        //cout << "mBTofPidTraitsIndex = " << track->bTofPidTraitsIndex() << endl;

        if(mCutManager->passTrackPhi(track,dca))     //shaowei
        {
            Float_t Mass2 = mCutManager->getMass2(track, pico); // shaowei
            //Float_t scale_nSigma_factor = ConstManager::mSigScaleMap[event->energy()];
            Float_t scale_nSigma_factor = 1.0;    // shaowei
            Float_t Polarity = static_cast<Float_t>(track->charge());
            Float_t momentum = track->pMom().Mag();
            Float_t Mass2_low;
            Float_t Mass2_up;
            /*if(momentum < 0.5)
            {
                Mass2_low = 0.4*0.4;
                Mass2_up = 0.6*0.6;
            }
            if(momentum >= 0.5)
            {
                Mass2_low = 0.277205 - 0.0812931*momentum;
                Mass2_up = 0.215517 + 0.076801*momentum;
            }
            if(momentum < 1.5)
            {
                Mass2_low = 0.4*0.4;
                Mass2_up = 0.6*0.6;
            }
            if(momentum >= 1.5)
            {
                Mass2_low = 0.15 * 0.15;
                Mass2_up = 0.6 * 0.6;
            }*/
            Mass2_low = 0.4*0.4;
            Mass2_up = 0.6*0.6;

            Int_t charge = 0; // k+
            if(Polarity < 0) charge = 1; // k-


            //if(mCutManager->passSigKaonCut(track, scale_nSigma_factor))
            //if(mCutManager->isKaon(pico, track))
            if(mCutManager->passSigKaonCut(track, scale_nSigma_factor))
            {
               if(
                        //((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0) // dE/dx + ToF
                        /*(momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
                        || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)*/
                            ((momentum <= 0.65 && Mass2 < -10) || (Mass2 > 0 && ((momentum < 1.5 && Mass2 > 0.16 && Mass2 < 0.36) || (momentum >= 1.5 && Mass2 > 0.125 && Mass2 < 0.36)) ))  &&
                            (
                             ((Mass2 < -10 && track->nSigmaKaon() < 2.5 && track->nSigmaKaon() > -(1.5)) || (Mass2 > 0.16 && Mass2 < 0.36)) 
                            )
                  )
                {
                    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
                    //mMass2[key].push_back(static_cast<Float_t>(mCutManager->getMass2(track, tof))); // mass2
                    mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
                    mDca[key].push_back(static_cast<Float_t>(dca*track->charge())); // dca*charge //shaowei
                    mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
                    mNSigmaKaon[key].push_back(static_cast<Float_t>((track->nSigmaKaon())*scale_nSigma_factor)); // nSigmaKaon
                    //mHelix_Kaon[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom 
                    mHelix_Kaon[key].push_back(static_cast<StPicoPhysicalHelix>(StPicoPhysicalHelix(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom 
                    mMomentum[key].push_back(static_cast<Float_t>(track->pMom().Mag()));// get helix from the pMom 
                    FillKaon(pico, track, mConfigs.y_mid);

               }
            }
	//cout << "mix event v0 test end" << endl;
        }
    }

    mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

    if(Flag_ME == 0) // same event
    {
	//cout << "do phi begin" << endl;
        doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
	//cout << "do phi end" << endl;
        clear_phi(cent9,Bin_vz,Bin_Psi2);
    }

    if(Flag_ME == 1) // mix event
    {
        if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == ConstManager::Buffer_depth)
        {
            doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
            clear_phi(cent9,Bin_vz,Bin_Psi2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------

// pass event information from Maker
void FlowV0::clearEvent()
{
    mNumber_prim = 0;
    mNumber_non_prim = 0;
    mNumber_Tof_match = 0;

    for(Int_t j = 0; j < ConstManager::EtaGap_total; j++)
    {
        mQVector1Full[j].Set(-999.9,-999.9);
        mQVector2East[j].Set(-999.9,-999.9);
        mQVector2West[j].Set(-999.9,-999.9);
        mQVector3East[j].Set(-999.9,-999.9);
        mQVector3West[j].Set(-999.9,-999.9);
    }
}

void FlowV0::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match)
{
    mNumber_prim = N_prim;
    mNumber_non_prim = N_non_prim;
    mNumber_Tof_match = N_Tof_match;
}

void FlowV0::passEventPlane1Full(TVector2 Q1Full_0, TVector2 Q1Full_1, TVector2 Q1Full_2, TVector2 Q1Full_3)
{
    mQVector1Full[0] = Q1Full_0;
    mQVector1Full[1] = Q1Full_1;
    mQVector1Full[2] = Q1Full_2;
    mQVector1Full[3] = Q1Full_3;
}

void FlowV0::passEventPlane2East(TVector2 Q2East_0, TVector2 Q2East_1, TVector2 Q2East_2, TVector2 Q2East_3)
{
    mQVector2East[0] = Q2East_0;
    mQVector2East[1] = Q2East_1;
    mQVector2East[2] = Q2East_2;
    mQVector2East[3] = Q2East_3;
}

void FlowV0::passEventPlane2West(TVector2 Q2West_0, TVector2 Q2West_1, TVector2 Q2West_2, TVector2 Q2West_3)
{
    mQVector2West[0] = Q2West_0;
    mQVector2West[1] = Q2West_1;
    mQVector2West[2] = Q2West_2;
    mQVector2West[3] = Q2West_3;
}

void FlowV0::passEventPlane3East(TVector2 Q3East_0, TVector2 Q3East_1, TVector2 Q3East_2, TVector2 Q3East_3)
{
    mQVector3East[0] = Q3East_0;
    mQVector3East[1] = Q3East_1;
    mQVector3East[2] = Q3East_2;
    mQVector3East[3] = Q3East_3;
}

void FlowV0::passEventPlane3West(TVector2 Q3West_0, TVector2 Q3West_1, TVector2 Q3West_2, TVector2 Q3West_3)
{
    mQVector3West[0] = Q3West_0;
    mQVector3West[1] = Q3West_1;
    mQVector3West[2] = Q3West_2;
    mQVector3West[3] = Q3West_3;
}

void FlowV0::passNumTrackEast(Int_t NumTrackEast_0, Int_t NumTrackEast_1, Int_t NumTrackEast_2, Int_t NumTrackEast_3)
{
    mTrackEtaEast[0] = NumTrackEast_0;
    mTrackEtaEast[1] = NumTrackEast_1;
    mTrackEtaEast[2] = NumTrackEast_2;
    mTrackEtaEast[3] = NumTrackEast_3;
}

void FlowV0::passNumTrackWest(Int_t NumTrackWest_0, Int_t NumTrackWest_1, Int_t NumTrackWest_2, Int_t NumTrackWest_3)
{
    mTrackEtaWest[0] = NumTrackWest_0;
    mTrackEtaWest[1] = NumTrackWest_1;
    mTrackEtaWest[2] = NumTrackWest_2;
    mTrackEtaWest[3] = NumTrackWest_3;
}
void FlowV0::FillKaon(StPicoDst *pico, StPicoTrack *track, Double_t y_mid)
{
  Double_t d_pT = track->pPt();
  Double_t d_px = track->pMom().X();
  Double_t d_py = track->pMom().Y();
  Double_t d_pz = track->pMom().Z();
  Double_t d_eta = track->pMom().Eta();
  Double_t d_phi = track->pMom().Phi();
  Double_t d_mom = track->pMom().Mag();
  //Double_t d_dEdx = track->dEdx();
  Short_t  s_charge = track->charge();
  TLorentzVector ltrack;
  ltrack.SetXYZM(d_px,d_py,d_pz,ConstManager::mMassKaon);
  Double_t mRapidity = ltrack.Rapidity();

  Double_t d_m2 = -999.0;
  Double_t d_tofBeta = -999.0;
  Int_t trackTofIndex = track->bTofPidTraitsIndex();
  if(trackTofIndex >= 0)
    d_tofBeta = pico->btofPidTraits(trackTofIndex)->btofBeta();
  if(d_tofBeta != -999.0)
    {
      d_m2 = d_mom*d_mom*( (1.0 / (d_tofBeta*d_tofBeta)) - 1.0 );
    }
  if(s_charge > 0)
  {
    h_eta_kp->Fill(d_eta);
    h_phi_kp->Fill(d_phi);
    h_pT_kp->Fill(d_pT);
    h_dndy_kp->Fill(mRapidity);
    h2_pT_vs_yCM_kp->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_kp->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_kp->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_kp->Fill(s_charge*d_mom, d_m2);
  }
  if(s_charge < 0)
  {
    h_eta_km->Fill(d_eta);
    h_phi_km->Fill(d_phi);
    h_pT_km->Fill(d_pT);
    h_dndy_km->Fill(mRapidity);
    h2_pT_vs_yCM_km->Fill(mRapidity - y_mid, d_pT);
    h2_dEdx_vs_qp_km->Fill(s_charge*d_mom, track->dEdx());
    h2_beta_vs_qp_km->Fill(s_charge*d_mom, 1.0/d_tofBeta);
    h2_m2_vs_qp_km->Fill(s_charge*d_mom, d_m2);
  }
}
//------------------------------------------------------------------------------------------------------------------
