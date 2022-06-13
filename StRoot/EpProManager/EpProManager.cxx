#include "EpProManager.h"
#include "../ConstManager/ConstManager.h"
#include "../Run/run.h"
//#include "StTriFlowConstants.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"

ClassImp(EpProManager)

    TString EpProManager::mVStr[2] = {"pos","neg"};

//---------------------------------------------------------------------------------

EpProManager::EpProManager()
{
}

//---------------------------------------------------------------------------------

EpProManager::~EpProManager()
{
}


//----------------------------------------------------------------------------

void EpProManager::InitEP()
{
    h2_RingIDvsEta_center = new TH2D("ringIDvsEta_center","rungIDvsEta_center",16,0,16,40,-6,-2);
    h2_RingIDvsEta_random = new TH2D("ringIDvsEta_random","rungIDvsEta_random",16,0,16,40,-6,-2);
    h_phi_center = new TH1D("h_phi_center","h_phi_center",630,-3.15,3.15);
    h_eta_center = new TH1D("h_eta_center","h_eta_center",300,-6,6);
    h_phi_random = new TH1D("h_phi_random","h_phi_random",630,-3.15,3.15);
    h_eta_random = new TH1D("h_eta_random","h_eta_random",300,-6,6);

    p_sector_eta=new TProfile2D("p_sector_eta", "", 300, -6.0, 6.0, 12, 0.5, 12.5);
    // recenter parameter for EPD EP
    /*for(int i=0; i<7; i++){
        p_EPDQx[i] = new TProfile2D(Form("EPDQx_ring%d_recen",i), "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
        p_EPDQy[i] = new TProfile2D(Form("EPDQy_ring%d_recen",i), "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    }
    // recenter parameter for TPC EP
    p_TPCqx_A = new TProfile2D("TPCqx_A_recen", "",  ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,500, -0.5, 499.5);
    p_TPCqy_A = new TProfile2D("TPCqy_A_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    p_TPCqx_B = new TProfile2D("TPCqx_B_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);
    p_TPCqy_B = new TProfile2D("TPCqy_B_recen", "", ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 500, -0.5, 499.5);*/
    v1EtaWt = new TH2D("v1EtaWt","v1EtaWt",500,1.5,6.5,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS);
    for (int ix=1; ix<501; ix++){
      for (int iy=1; iy<10; iy++){
        double eta = v1EtaWt->GetXaxis()->GetBinCenter(ix);
        v1EtaWt->SetBinContent(ix,iy,_lin[iy-1]*eta+_cub[iy-1]*pow(eta,3)+_fif[iy-1]*pow(eta,5));
      }
    }
    

    for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
    {
    	//
    	float param1 = _subEventParams[sub][0];// - COMrapidity;
    	float param2 = _subEventParams[sub][1];// - COMrapidity;
    	float param1_tpc = _subEventParams_tpc[sub][0];// - COMrapidity;
    	float param2_tpc = _subEventParams_tpc[sub][1];// - COMrapidity;
    	
    	TString form = Form(" - #sub event #in %c (%1.1f, %1.1f)",_subEventModes[sub],param1 ,param2);
    	TString form_tpc = Form(" - #sub event #in %c (%1.1f, %1.1f)",_subEventModes[sub],param1_tpc ,param2_tpc);
    	
    	h_psi1_epd_ABCD_raw_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_raw_sub_%d",sub + 1), "#Psi_{1}^{EPD} distribution (raw)" + form,1024,-7.0,7.0);
    	h_psi1_epd_ABCD_raw_wt_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_raw_wt_sub_%d",sub + 1), "#Psi_{1}^{EPD} w/ eta weighting distribution (raw)" + form,1024,-7.0,7.0);
    	h_psi2_tpc_AB_raw_sub[sub] = new TH1F(Form("h_psi2_tpc_AB_raw_sub_%d",sub + 1), "#Psi_{2}^{TPC} distribution (raw)" + form_tpc,1024,-7.0,7.0);
    	h_psi1_epd_ABCD_recentered_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_recentered_sub_%d",sub + 1), "#Psi_{1}^{EPD} distribution (recentered)" + form,1024,-7.0,7.0);
    	h_psi1_epd_ABCD_recentered_wt_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_recentered_wt_sub_%d",sub + 1), "#Psi_{1}^{EPD} w/ eta weighting distribution (recentered)" + form,1024,-7.0,7.0);
    	h_psi2_tpc_AB_recentered_sub[sub] = new TH1F(Form("h_psi2_tpc_AB_recentered_sub_%d",sub + 1), "#Psi_{2}^{TPC} distribution (recentered)" + form_tpc,1024,-7.0,7.0);
    	h_psi1_epd_ABCD_shifted_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_shifted_sub_%d",sub + 1), "#Psi_{1}^{EPD} distribution (shifted)" + form,1024,-7.0,7.0);
    	h_psi1_epd_ABCD_shifted_wt_sub[sub] = new TH1F(Form("h_psi1_epd_ABCD_shifted_wt_sub_%d",sub + 1), "#Psi_{1}^{EPD} w/ eta weighting distribution (shifted)" + form,1024,-7.0,7.0);
    	h_psi2_tpc_AB_shifted_sub[sub] = new TH1F(Form("h_psi2_tpc_AB_shifted_sub_%d",sub + 1), "#Psi_{2}^{TPC} distribution (shifted)" + form_tpc,1024,-7.0,7.0);
    	
    	hQxSub[sub] = new TH1F(Form("hQxSub_{%d}",sub + 1), "Q_{x} distribution",10000,-100,100);	
    	hQySub[sub] = new TH1F(Form("hQySub_{%d}",sub + 1), "Q_{x} distribution",10000,-100,100);	
    	hQVectorSub[sub] = new TH2F(Form("hqVectorSub_{%d}",sub + 1),"#vec{Q} distribution",200,-100,100,200,-100,100);
    	hQxSubRec[sub] = new TH1F(Form("hQxSubRec_%d",sub + 1), "Q_{x}Rec distribution",1000,-100,100);	
    	hQySubRec[sub] = new TH1F(Form("hQySubRec_%d",sub + 1), "Q_{x}Rec distribution",1000,-100,100);	
    	hQVectorSubRec[sub] = new TH2F(Form("hqVectorSubRec_%d",sub + 1),"#vec{Q}Rec distribution",20,-100,100,20,-100,100);
    	p_mq1x_epd_ABCD_sub[sub] = new TProfile2D(Form("p_mq1x_epd_ABCD_sub_%d",sub + 1),"mQx1 for recenter correction Qx" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    	p_mq1y_epd_ABCD_sub[sub] = new TProfile2D(Form("p_mq1y_epd_ABCD_sub_%d",sub + 1),"mQy1 for recenter correction Qy" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    	p_mq1x_epd_ABCD_wt_sub[sub] = new TProfile2D(Form("p_mq1x_epd_ABCD_wt_sub_%d",sub + 1),"mQx1 for recenter correction Qx EPD w/ wt" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    	p_mq1y_epd_ABCD_wt_sub[sub] = new TProfile2D(Form("p_mq1y_epd_ABCD_wt_sub_%d",sub + 1),"mQy1 for recenter correction Qx EPD w/ wt" + form,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    	p_mq2x_tpc_AB_sub[sub] = new TProfile2D(Form("p_mq2x_tpc_AB_sub_%d",sub + 1),"mQx2 for recenter correction Qx" + form_tpc,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    	p_mq2y_tpc_AB_sub[sub] = new TProfile2D(Form("p_mq2y_tpc_AB_sub_%d",sub + 1),"mQy2 for recenter correction Qy" + form_tpc,ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
	// shift parameter
        p_shiftpar_sin_epd_ABCD_sub[sub] =new TProfile3D(Form("EPshiftpar_epd_ABCD_sub_%d_sin", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
        p_shiftpar_cos_epd_ABCD_sub[sub] =new TProfile3D(Form("EPshiftpar_epd_ABCD_sub_%d_cos", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
        p_shiftpar_sin_epd_ABCD_wt_sub[sub] =new TProfile3D(Form("EPshiftpar_epd_ABCD_wt_sub_%d_sin", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
        p_shiftpar_cos_epd_ABCD_wt_sub[sub] =new TProfile3D(Form("EPshiftpar_epd_ABCD_wt_sub_%d_cos", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
        p_shiftpar_sin_tpc_AB_sub[sub] =new TProfile3D(Form("EPshiftpar_tpc_AB_sub_%d_sin", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
        p_shiftpar_cos_tpc_AB_sub[sub] =new TProfile3D(Form("EPshiftpar_tpc_AB_sub_%d_cos", sub + 1),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    } // event plane Psi histograms
    
    h_psi1_epd_ABCD_raw_full = new TH1F("h_psi1_epd_ABCD_raw_full", "#Psi_{1}^{EPD-Full} distribution (raw)" ,1024,-7.0,7.0);
    h_psi1_epd_ABCD_raw_wt_full = new TH1F("h_psi1_epd_ABCD_raw_wt_full", "#Psi_{1}^{EPD-Full} w/ eta weighting distribution (raw)" ,1024,-7.0,7.0);
    h_psi2_tpc_AB_raw_full = new TH1F("h_psi2_tpc_AB_raw_full", "#Psi_{2}^{TPC-Full} distribution (raw)" ,1024,-7.0,7.0);
    h_psi1_epd_ABCD_raw_subs = new TH2F(Form("h_psi1_epd_ABCD_raw_subs_%d", 1), "#Psi distribution (raw) #Psi_{1}^{East}  vs. #Psi_{1}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi1_epd_ABCD_raw_wt_subs = new TH2F(Form("h_psi1_epd_ABCD_raw_wt_subs_%d", 1), "#Psi distribution (raw) w/ eta weighting #Psi_{1}^{EPD-East}  vs. #Psi_{1}^{EPD-West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi2_tpc_AB_raw_subs = new TH2F(Form("h_psi2_tpc_AB_raw_subs_%d", 1), "#Psi distribution (raw) #Psi_{2}^{TPC-East}  vs. #Psi_{2}^{TPC-West} " ,35,-3.5,3.5,35,-3.5,3.5);
    	
    h_psi1_epd_ABCD_recentered_full = new TH1F(Form("h_psi1_epd_ABCD_recentered_full"), "#Psi_{1}^{EPD-full} distribution (recentered)",1024,-7.0,7.0);
    h_psi1_epd_ABCD_recentered_subs = new TH2F(Form("h_psi1_epd_ABCD_recentered_subs"), "#Psi_{1} distribution (recentered) #Psi_{1}^{East}  vs. #Psi_{1}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi1_epd_ABCD_recentered_wt_full = new TH1F(Form("h_psi1_epd_ABCD_recentered_wt_full"), "#Psi_{1}^{EPD-full} w/ eta weighting distribution (recentered)",1024,-7.0,7.0);
    h_psi1_epd_ABCD_recentered_wt_subs = new TH2F(Form("h_psi1_epd_ABCD_recentered_wt_subs"), "#Psi_{1}^{EPD-full} w/ eat weighting distribution (recentered) #Psi_{1}^{East}  vs. #Psi_{1}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi2_tpc_AB_recentered_full = new TH1F(Form("h_psi2_tpc_AB_recentered_full"), "#Psi_{2}^{TPC-full} distribution (recentered)",1024,-7.0,7.0);
    h_psi2_tpc_AB_recentered_subs = new TH2F(Form("h_psi2_tpc_AB_recentered_subs"), "#Psi_{2} distribution (recentered) #Psi_{2}^{East}  vs. #Psi_{2}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi1_epd_ABCD_shifted_full = new TH1F(Form("h_psi1_epd_ABCD_shifted_full"), "#Psi_{1}^{EPD-full} distribution (shifted)",1024,-7.0,7.0);
    h_psi1_epd_ABCD_shifted_subs = new TH2F(Form("h_psi1_epd_ABCD_shifted_subs"), "#Psi_{1} distribution (shifted) #Psi_{1}^{East}  vs. #Psi_{1}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi1_epd_ABCD_shifted_wt_full = new TH1F(Form("h_psi1_epd_ABCD_shifted_wt_full"), "#Psi_{1}^{EPD-full} w/ eta weighting distribution (shifted)",1024,-7.0,7.0);
    h_psi1_epd_ABCD_shifted_wt_subs = new TH2F(Form("h_psi1_epd_ABCD_shifted_wt_subs"), "#Psi_{1}^{EPD-full} w/ eat weighting distribution (shifted) #Psi_{1}^{East}  vs. #Psi_{1}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    h_psi2_tpc_AB_shifted_full = new TH1F(Form("h_psi2_tpc_AB_shifted_full"), "#Psi_{2}^{TPC-full} distribution (shifted)",1024,-7.0,7.0);
    h_psi2_tpc_AB_shifted_subs = new TH2F(Form("h_psi2_tpc_AB_shifted_subs"), "#Psi_{2} distribution (shifted) #Psi_{2}^{East}  vs. #Psi_{2}^{West} " ,35,-3.5,3.5,35,-3.5,3.5);
    // parameter for recentering for full event planes 
    p_mq1x_epd_ABCD_full = new TProfile2D(Form("p_mq1x_epd_ABCD_full"),"mQx1 for recenter correction full Qx",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    p_mq1y_epd_ABCD_full = new TProfile2D(Form("p_mq1y_epd_ABCD_full"),"mQy1 for recenter correction Qy",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    p_mq1x_epd_ABCD_wt_full = new TProfile2D(Form("p_mq1x_epd_ABCD_wt_full"),"mQx1 w/ eta weight for recenter correction full Qx",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    p_mq1y_epd_ABCD_wt_full = new TProfile2D(Form("p_mq1y_epd_ABCD_wt_full"),"mQy1 w/ eta weight for recenter correction Qy",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    p_mq2x_tpc_AB_full = new TProfile2D(Form("p_mq2x_tpc_AB_full"),"mQx2 tpc for recenter correction full Qx",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    p_mq2y_tpc_AB_full = new TProfile2D(Form("p_mq2y_tpc_AB_full"),"mQy2 tpc  for recenter correction Qy",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, nrun, -0.5, -0.5+nrun);
    // shift parameter
    p_shiftpar_sin_epd_ABCD_full =new TProfile3D(Form("EPshiftpar_epd_ABCD_full_sin"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    p_shiftpar_cos_epd_ABCD_full =new TProfile3D(Form("EPshiftpar_epd_ABCD_full_cos"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    p_shiftpar_sin_epd_ABCD_wt_full =new TProfile3D(Form("EPshiftpar_epd_ABCD_wt_full_sin"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    p_shiftpar_cos_epd_ABCD_wt_full =new TProfile3D(Form("EPshiftpar_epd_ABCD_wt_full_cos"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    p_shiftpar_sin_tpc_AB_full =new TProfile3D(Form("EPshiftpar_tpc_AB_full_sin"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    p_shiftpar_cos_tpc_AB_full =new TProfile3D(Form("EPshiftpar_tpc_AB_full_cos"),"",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS, 20, -0.5, 19.5, nrun, -0.5, -0.5+nrun);
    // resolution plots
    p_r1_epd_ABCD_sub0_1 = new TProfile("p_r1_epd_ABCD_sub0_1","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    p_r1_epd_ABCD_sub0_1_sin = new TProfile("p_r1_epd_ABCD_sub0_1_sin","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    p_r1_epd_ABCD_wt_sub0_1 = new TProfile("p_r1_epd_ABCD_wt_sub0_1","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    p_r1_epd_ABCD_wt_sub0_1_sin = new TProfile("p_r1_epd_ABCD_wt_sub0_1_sin","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    p_r2_tpc_AB_sub0_1 = new TProfile("p_r2_tpc_AB_sub0_1","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    p_r2_tpc_AB_sub0_1_sin = new TProfile("p_r2_tpc_AB_sub0_1_sin","",ConstManager::CENT_BINS, ConstManager::FIRST_CENT, ConstManager::FIRST_CENT + ConstManager::CENT_BINS,-1.0,1.0,"");
    cout << "End of event plane Histograms" << endl;
}


/*void EpProManager::InitEP()
{
    //===============================                          
    //  Define Histograms                           
    //===============================                                            

    h2_RingIDvsEta_center = new TH2D("ringIDvsEta_center","rungIDvsEta_center",16,0,16,40,-6,-2);
    h2_RingIDvsEta_random = new TH2D("ringIDvsEta_random","rungIDvsEta_random",16,0,16,40,-6,-2);
    h_phi_center = new TH1D("h_phi_center","h_phi_center",630,-3.15,3.15);
    h_eta_center = new TH1D("h_eta_center","h_eta_center",300,-6,6);
    h_phi_random = new TH1D("h_phi_random","h_phi_random",630,-3.15,3.15);
    h_eta_random = new TH1D("h_eta_random","h_eta_random",300,-6,6);


    for(int i=0; i<9; i++){
        h_eta[i] = new TH1F(Form("eta_%d",i),Form("eta_%d",i),300,-2.5,0.5);
        h2_eta_pt[i] = new TH2D(Form("eta_pt_%d",i),Form("eta_pt_%d",i),300,-2.5,0.5,200,0.0,5.0);
    }


    cout << "End of vent plane Histograms" << endl;
 
    for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                TString ProName;
                // Event Plane method
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin2_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin2_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_mcos3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[i].Data(),j,k);
                p_msin3_East_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_mcos3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
                ProName = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[i].Data(),j,k);
                p_msin3_West_EP[i][j][k] = new TProfile2D(ProName.Data(),ProName.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

                // Scalor Product method
            }
        }

        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            TString ProName_Full;
            // Event Plane method
            ProName_Full = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName_Full = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin2_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

            ProName_Full = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_mcos3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality
            ProName_Full = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[i].Data(),k);
            p_msin3_Full_EP[i][k] = new TProfile2D(ProName_Full.Data(),ProName_Full.Data(),1600,-0.5,1599.5,9,-0.5,8.5); // x axis is RunIndex, y axis is Centrality

        }
    }

}
*/

//----------------------------------------------------------------------------
// Event Plane method
void EpProManager::FillEpdQa(Int_t iring, Double_t phi_epd_center, Double_t eta_epd_center,Double_t phi_epd_random, Double_t eta_epd_random)
{
  h_phi_center -> Fill(phi_epd_center);
  h_eta_center -> Fill(eta_epd_center);
  h_phi_random -> Fill(phi_epd_random);
  h_eta_random -> Fill(eta_epd_random);
  h2_RingIDvsEta_center -> Fill(iring,eta_epd_center);
  h2_RingIDvsEta_random -> Fill(iring,eta_epd_random);

}
void EpProManager::FillEpdMip(Double_t eta_epd_center, Int_t position, Double_t mip)
{
  p_sector_eta->Fill(eta_epd_center, position, mip);
}
void EpProManager::FillEpdQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_EPDQx[isub]->Fill(centnumber, runindex, qx);
  p_EPDQy[isub]->Fill(centnumber, runindex, qy);
}

void EpProManager::FillTpcAQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track)
{
  Double_t d_pT = track->pPt();
  Double_t d_phi = track->pMom().Phi();
  Double_t qx_TPC_A = d_pT*TMath::Cos(1.0*d_phi);
  Double_t qy_TPC_A = d_pT*TMath::Sin(1.0*d_phi);
  p_TPCqx_A->Fill(centnumber, runindex, qx_TPC_A);
  p_TPCqy_A->Fill(centnumber, runindex, qy_TPC_A);

}
void EpProManager::FillTpcBQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track)
{
  Double_t d_pT = track->pPt();
  Double_t d_phi = track->pMom().Phi();
  Double_t qx_TPC_A = d_pT*TMath::Cos(1.0*d_phi);
  Double_t qy_TPC_A = d_pT*TMath::Sin(1.0*d_phi);
  p_TPCqx_B->Fill(centnumber, runindex, qx_TPC_A);
  p_TPCqy_B->Fill(centnumber, runindex, qy_TPC_A);

}
void EpProManager::FillSubEpQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq1x_epd_ABCD_sub[isub]->Fill(centnumber, runindex, qx);
  p_mq1y_epd_ABCD_sub[isub]->Fill(centnumber, runindex, qy);
  hQxSub[isub]->Fill(qx);
  hQySub[isub]->Fill(qy);
  hQVectorSub[isub]->Fill(qx,qy);
}

void EpProManager::FillSubEpQvec_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq1x_epd_ABCD_full->Fill(centnumber, runindex, qx);
  p_mq1y_epd_ABCD_full->Fill(centnumber, runindex, qy);
}

void EpProManager::FillSubEpQvec_wt(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq1x_epd_ABCD_wt_sub[isub]->Fill(centnumber, runindex, qx);
  p_mq1y_epd_ABCD_wt_sub[isub]->Fill(centnumber, runindex, qy);
}
    
void EpProManager::FillSubEpQvec_wt_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq1x_epd_ABCD_wt_full->Fill(centnumber, runindex, qx);
  p_mq1y_epd_ABCD_wt_full->Fill(centnumber, runindex, qy);
}

void EpProManager::FillSubEpQvec_tpc(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq2x_tpc_AB_sub[isub]->Fill(centnumber, runindex, qx);
  p_mq2y_tpc_AB_sub[isub]->Fill(centnumber, runindex, qy);
}

void EpProManager::FillSubEpQvec_tpc_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  p_mq2x_tpc_AB_full->Fill(centnumber, runindex, qx);
  p_mq2y_tpc_AB_full->Fill(centnumber, runindex, qy);
}

void EpProManager::FillPsiRaw(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_raw_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRaw_wt(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_raw_wt_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRawFull( Double_t psi)
{
  h_psi1_epd_ABCD_raw_full->Fill(psi);
}

void EpProManager::FillPsiRawFull_wt( Double_t psi)
{
  h_psi1_epd_ABCD_raw_wt_full->Fill(psi);
}

void EpProManager::FillPsiRaw_tpc(Int_t isub, Double_t psi)
{
  h_psi2_tpc_AB_raw_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRawFull_tpc( Double_t psi)
{
  h_psi2_tpc_AB_raw_full->Fill(psi);
}
void EpProManager::FillPsiRawSubs(Double_t psi_east, Double_t psi_west)
{
  h_psi1_epd_ABCD_raw_subs->Fill(psi_east, psi_west);
}

void EpProManager::FillPsiRawSubs_wt(Double_t psi_east, Double_t psi_west)
{
  h_psi1_epd_ABCD_raw_wt_subs->Fill(psi_east, psi_west);
}
void EpProManager::FillPsiRawSubs_tpc(Double_t psi_east, Double_t psi_west)
{
  h_psi2_tpc_AB_raw_subs->Fill(psi_east, psi_west);
}

// Fill recenter histograms

void EpProManager::FillSubEpQvecRec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy)
{
  hQxSubRec[isub]->Fill(qx);
  hQySubRec[isub]->Fill(qy);
  hQVectorSubRec[isub]->Fill(qx,qy);
}

void EpProManager::FillPsiRec(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_recentered_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRec_wt(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_recentered_wt_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRecFull_wt( Double_t psi)
{
  h_psi1_epd_ABCD_recentered_wt_full->Fill(psi);
}

void EpProManager::FillPsiRecFull_tpc( Double_t psi)
{
  h_psi2_tpc_AB_recentered_full->Fill(psi);
}
void EpProManager::FillPsiRec_tpc(Int_t isub, Double_t psi)
{
  h_psi2_tpc_AB_recentered_sub[isub]->Fill(psi);
}

void EpProManager::FillPsiRecFull( Double_t psi)
{
  h_psi1_epd_ABCD_recentered_full->Fill(psi);
}

void EpProManager::FillPsiRecSubs(Double_t psi_east, Double_t psi_west)
{
  h_psi1_epd_ABCD_recentered_subs->Fill(psi_east, psi_west);
}
void EpProManager::FillPsiRecSubs_wt(Double_t psi_east, Double_t psi_west)
{
  h_psi1_epd_ABCD_recentered_wt_subs->Fill(psi_east, psi_west);
}
void EpProManager::FillPsiRecSubs_tpc(Double_t psi_east, Double_t psi_west)
{
  h_psi2_tpc_AB_recentered_subs->Fill(psi_east, psi_west);
}

// Fill shift histograms
void EpProManager::FillPsiShift(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_shifted_sub[isub]->Fill(psi);
}
void EpProManager::FillPsiShiftFull( Double_t psi)
{
  h_psi1_epd_ABCD_shifted_full->Fill(psi);
}
void EpProManager::FillPsiShiftSubs( Double_t psi0, Double_t psi1)
{
  h_psi1_epd_ABCD_shifted_subs->Fill(psi0, psi1);
}

void EpProManager::FillPsiShift_wt(Int_t isub, Double_t psi)
{
  h_psi1_epd_ABCD_shifted_wt_sub[isub]->Fill(psi);
}
void EpProManager::FillPsiShiftFull_wt( Double_t psi)
{
  h_psi1_epd_ABCD_shifted_wt_full->Fill(psi);
}
void EpProManager::FillPsiShiftSubs_wt( Double_t psi0, Double_t psi1)
{
  h_psi1_epd_ABCD_shifted_wt_subs->Fill(psi0, psi1);
}

void EpProManager::FillPsiShift_tpc(Int_t isub, Double_t psi)
{
  h_psi2_tpc_AB_shifted_sub[isub]->Fill(psi);
}
void EpProManager::FillPsiShiftFull_tpc( Double_t psi)
{
  h_psi2_tpc_AB_shifted_full->Fill(psi);
}
void EpProManager::FillPsiShiftSubs_tpc( Double_t psi0, Double_t psi1)
{
  h_psi2_tpc_AB_shifted_subs->Fill(psi0, psi1);
}
// Fill shiftpar for the shift correction
void EpProManager::FillSubEpShiftpar(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_epd_ABCD_sub[isub]->Fill(centnumber, iorder, runindex, sin(1.0*(iorder+1)*psi));
  p_shiftpar_cos_epd_ABCD_sub[isub]->Fill(centnumber, iorder, runindex, cos(1.0*(iorder+1)*psi));
}
void EpProManager::FillFullEpShiftpar(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_epd_ABCD_full->Fill(centnumber, iorder, runindex, sin(1.0*(iorder+1)*psi));
  p_shiftpar_cos_epd_ABCD_full->Fill(centnumber, iorder, runindex, cos(1.0*(iorder+1)*psi));
}

void EpProManager::FillSubEpShiftpar_wt(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_epd_ABCD_wt_sub[isub]->Fill(centnumber, iorder, runindex, sin(1.0*(iorder+1)*psi));
  p_shiftpar_cos_epd_ABCD_wt_sub[isub]->Fill(centnumber, iorder, runindex, cos(1.0*(iorder+1)*psi));
}
void EpProManager::FillFullEpShiftpar_wt(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_epd_ABCD_wt_full->Fill(centnumber, iorder, runindex, sin(1.0*(iorder+1)*psi));
  p_shiftpar_cos_epd_ABCD_wt_full->Fill(centnumber, iorder, runindex, cos(1.0*(iorder+1)*psi));
}

void EpProManager::FillSubEpShiftpar_tpc(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_tpc_AB_sub[isub]->Fill(centnumber, iorder, runindex, sin(2.0*(iorder+1)*psi));
  p_shiftpar_cos_tpc_AB_sub[isub]->Fill(centnumber, iorder, runindex, cos(2.0*(iorder+1)*psi));
}
void EpProManager::FillFullEpShiftpar_tpc(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi)
{
  p_shiftpar_sin_tpc_AB_full->Fill(centnumber, iorder, runindex, sin(2.0*(iorder+1)*psi));
  p_shiftpar_cos_tpc_AB_full->Fill(centnumber, iorder, runindex, cos(2.0*(iorder+1)*psi));
}

void EpProManager::FillPsi1ResolutionEpd(Int_t centnumber, Double_t psi0, Double_t psi1)
{
    double r1_Sub0_1 = TMath::Cos(1.0*(psi0 - psi1));
    double r1_Sub0_1_sin = TMath::Sin(1.0*(psi0 - psi1));

    p_r1_epd_ABCD_sub0_1 -> Fill(centnumber, r1_Sub0_1);
    p_r1_epd_ABCD_sub0_1_sin -> Fill(centnumber, r1_Sub0_1_sin);

}
void EpProManager::FillPsi1ResolutionEpdWt(Int_t centnumber, Double_t psi0, Double_t psi1)
{
    double r1_Sub0_1 = TMath::Cos(1.0*(psi0 - psi1));
    double r1_Sub0_1_sin = TMath::Sin(1.0*(psi0 - psi1));

    p_r1_epd_ABCD_wt_sub0_1 -> Fill(centnumber, r1_Sub0_1);
    p_r1_epd_ABCD_wt_sub0_1_sin -> Fill(centnumber, r1_Sub0_1_sin);

}
void EpProManager::FillPsi2ResolutionTpc(Int_t centnumber, Double_t psi0, Double_t psi1)
{
    double r2_Sub0_1 = TMath::Cos(2.0*(psi0 - psi1));
    double r2_Sub0_1_sin = TMath::Sin(2.0*(psi0 - psi1));

    p_r2_tpc_AB_sub0_1 -> Fill(centnumber, r2_Sub0_1);
    p_r2_tpc_AB_sub0_1_sin -> Fill(centnumber, r2_Sub0_1_sin);

}
/*void EpProManager::FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_East_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

/*void EpProManager::FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k) // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_West_EP[i][j][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

/*void EpProManager::FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k) // i = vertex pos/neg, k = ShiftOrder
{
    const Float_t cos2 = Psi2Vector.X();
    const Float_t sin2 = Psi2Vector.Y();
    const Float_t cos3 = Psi3Vector.X();
    const Float_t sin3 = Psi3Vector.Y();

    p_mcos2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos2);
    p_msin2_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin2);

    p_mcos3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)cos3);
    p_msin3_Full_EP[i][k]->Fill((Double_t)RunIndex,(Double_t)Cent9,(Double_t)sin3);
}*/

//----------------------------------------------------------------------------
Double_t EpProManager::GetEtaWeight(Int_t centrality, Double_t eta)
{
	int v1etaBin = (int)v1EtaWt->GetXaxis()->FindBin(abs(eta));
	int centBin = (int)v1EtaWt->GetYaxis()->FindBin(centrality);
        double v1EtaWeight = (double)v1EtaWt->GetBinContent(v1etaBin,centBin);
 	return v1EtaWeight;
}
//----------------------------------------------------------------------------

void EpProManager::WriteEP()
{
  h_phi_center -> Write();
  h_eta_center -> Write();
  h_phi_random -> Write();
  h_eta_random -> Write();
  h2_RingIDvsEta_center -> Write();
  h2_RingIDvsEta_random -> Write();
  p_sector_eta -> Write();
  v1EtaWt->Write();
  /*for(int i=0; i<7; i++){
      p_EPDQx[i]->Write();
      p_EPDQy[i]->Write();
  }
  p_TPCqx_A->Write();
  p_TPCqy_A->Write();
  p_TPCqx_B->Write();
  p_TPCqy_B->Write();*/
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    hQxSub[sub]->Write();
    hQySub[sub]->Write();
    hQVectorSub[sub]->Write();
  }
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_raw_sub[sub] -> Write();
  }
  h_psi1_epd_ABCD_raw_full->Write();
  h_psi1_epd_ABCD_raw_subs->Write();
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_recentered_sub[sub] -> Write();
  }
  h_psi1_epd_ABCD_recentered_full->Write();
  h_psi1_epd_ABCD_recentered_subs->Write();
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_shifted_sub[sub] -> Write();
  }
  h_psi1_epd_ABCD_shifted_full->Write();
  h_psi1_epd_ABCD_shifted_subs->Write();
  
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_raw_wt_sub[sub] ->Write();
  }
  h_psi1_epd_ABCD_raw_wt_full->Write();
  h_psi1_epd_ABCD_raw_wt_subs->Write();

  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_recentered_wt_sub[sub] -> Write();
  }
  h_psi1_epd_ABCD_recentered_wt_full->Write();
  h_psi1_epd_ABCD_recentered_wt_subs->Write();

  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi1_epd_ABCD_shifted_wt_sub[sub] -> Write();
  }
  h_psi1_epd_ABCD_shifted_wt_full->Write();
  h_psi1_epd_ABCD_shifted_wt_subs->Write();

  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi2_tpc_AB_raw_sub[sub] ->Write();
  }
  h_psi2_tpc_AB_raw_full->Write();
  h_psi2_tpc_AB_raw_subs->Write();
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi2_tpc_AB_recentered_sub[sub] -> Write();
  }
  h_psi2_tpc_AB_recentered_full->Write();
  h_psi2_tpc_AB_recentered_subs->Write();
  
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    h_psi2_tpc_AB_shifted_sub[sub] -> Write();
  }
  h_psi2_tpc_AB_shifted_full->Write();
  h_psi2_tpc_AB_shifted_subs->Write();
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_mq1x_epd_ABCD_sub[sub]->Write();
    p_mq1y_epd_ABCD_sub[sub]->Write();
  }
    p_mq1x_epd_ABCD_full->Write();
    p_mq1y_epd_ABCD_full->Write();
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_mq1x_epd_ABCD_wt_sub[sub]->Write();
    p_mq1y_epd_ABCD_wt_sub[sub]->Write();
  }
    p_mq1x_epd_ABCD_wt_full->Write();
    p_mq1y_epd_ABCD_wt_full->Write();
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_mq2x_tpc_AB_sub[sub]->Write();
    p_mq2y_tpc_AB_sub[sub]->Write();
  }
    p_mq2x_tpc_AB_full->Write();
    p_mq2y_tpc_AB_full->Write();
 // shift parameters 
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_shiftpar_sin_epd_ABCD_sub[sub]->Write();
    p_shiftpar_cos_epd_ABCD_sub[sub]->Write();
  }
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_shiftpar_sin_epd_ABCD_wt_sub[sub]->Write();
    p_shiftpar_cos_epd_ABCD_wt_sub[sub]->Write();
  }
  for (int sub = 0; sub < _numSubEvents; sub++) // event plane Psi histograms
  {
    p_shiftpar_sin_tpc_AB_sub[sub]->Write();
    p_shiftpar_cos_tpc_AB_sub[sub]->Write();
  }
    p_shiftpar_sin_epd_ABCD_full->Write();
    p_shiftpar_cos_epd_ABCD_full->Write();
    p_shiftpar_sin_epd_ABCD_wt_full->Write();
    p_shiftpar_cos_epd_ABCD_wt_full->Write();
    p_shiftpar_sin_tpc_AB_full->Write();
    p_shiftpar_cos_tpc_AB_full->Write();
    p_r1_epd_ABCD_sub0_1->Write();
    p_r1_epd_ABCD_sub0_1_sin->Write();
    p_r1_epd_ABCD_wt_sub0_1->Write();
    p_r1_epd_ABCD_wt_sub0_1_sin->Write();
    p_r2_tpc_AB_sub0_1->Write();
    p_r2_tpc_AB_sub0_1_sin->Write();
    /*for(Int_t i = 0; i < 2; i++) // vertex pos/neg
    {
        for(Int_t j = 0; j < 4; j++) // eta_gap
        {
            for(Int_t k = 0; k < 5; k++) // Shift Order
            {
                // Event Plane method
                p_mcos2_East_EP[i][j][k]->Write();
                p_msin2_East_EP[i][j][k]->Write();
                p_mcos2_West_EP[i][j][k]->Write();
                p_msin2_West_EP[i][j][k]->Write();

                p_mcos3_East_EP[i][j][k]->Write();
                p_msin3_East_EP[i][j][k]->Write();
                p_mcos3_West_EP[i][j][k]->Write();
                p_msin3_West_EP[i][j][k]->Write();

            }
        }
        for(Int_t k = 0; k < 5; k++) // Shift Order
        {
            // Event Plane method
            p_mcos2_Full_EP[i][k]->Write();
            p_msin2_Full_EP[i][k]->Write();
            p_mcos3_Full_EP[i][k]->Write();
            p_msin3_Full_EP[i][k]->Write();

        }
    }*/
}

//----------------------------------------------------------------------------


