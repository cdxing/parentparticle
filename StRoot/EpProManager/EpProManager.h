#ifndef EpProManager_h
#define EpProManager_h

#include "TVector2.h"
#include "TString.h"
#include "StMessMgr.h"
#include "../StPicoEvent/StPicoTrack.h"

class TProfile3D;
class TProfile2D;
class TProfile;
class TH1F;
class TH1D;
class TH2F;
class TH2D;

// sub event plane configurations
const int _numSubEvents = 2;
const char _subEventModes[_numSubEvents] = {'e','e'};
const float _subEventParams[_numSubEvents][2] =
{{-5.4,-1.6}, //EPD-East
 {1.6,5.4}, //EPD-West
 }; 
const float _subEventParams_tpc[_numSubEvents][2] =
{{-1.5,-0.05}, //TPC-East
 {0.05,1.5}, //TPC-West
 }; 
// eta weight from Zuowen's study:
// https://drupal.star.bnl.gov/STAR/system/files/Directed_Flow_19p6GeV_0.pdf
const  double _lin[9] = {-2.23911, -2.57653, -2.51053, -2.32201, -1.95298, -1.46519, -0.91108, -0.41001, 0.035437};
const  double _cub[9] = {0.364133, 0.339185, 0.311601, 0.291423, 0.24892, 0.197276, 0.128563, 0.069993, 0.013928};
const  double _fif[9] = {-0.01212,  -0.00716,  -0.00465,  -0.0034, -0.00189, -0.00086,  0.000363,  0.000486,  0.00046};


class EpProManager
{
  public:
    EpProManager();
    ~EpProManager();


    // Recenter and  Shift Correction
    void InitEP();
    void WriteEP();
    // Event Plane method
    // Fill plots
    void FillEpdQa(Int_t iring, Double_t phi_epd_center, Double_t eta_epd_center,Double_t phi_epd_random, Double_t eta_epd_random); 
    void FillEpdMip(Double_t eta_epd_center, Int_t position, Double_t mip); 
    void FillEpdQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy); 
    void FillTpcAQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track);
    void FillTpcBQvec(Int_t centnumber, Int_t runindex, StPicoTrack *track);
    void FillSubEpQvec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    void FillSubEpQvec_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    
    void FillSubEpQvec_wt(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    void FillSubEpQvec_wt_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    
    void FillSubEpQvec_tpc(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    void FillSubEpQvec_tpc_full(Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    
    void FillPsiRaw(Int_t isub, Double_t psi); 
    void FillPsiRaw_wt(Int_t isub, Double_t psi); 
    void FillPsiRaw_tpc(Int_t isub, Double_t psi); 
    void FillPsiRawFull( Double_t psi);
    void FillPsiRawFull_wt( Double_t psi);
    void FillPsiRawFull_tpc( Double_t psi);
    void FillPsiRawSubs(Double_t psi_east, Double_t psi_west);
    void FillPsiRawSubs_wt(Double_t psi_east, Double_t psi_west);
    void FillPsiRawSubs_tpc(Double_t psi_east, Double_t psi_west);
    void FillPsiRec(Int_t isub, Double_t psi); 
    void FillPsiRec_wt(Int_t isub, Double_t psi); 
    void FillPsiRec_tpc(Int_t isub, Double_t psi); 
    void FillPsiRecFull( Double_t psi);
    void FillPsiRecFull_wt( Double_t psi);
    void FillPsiRecFull_tpc( Double_t psi);
    void FillPsiRecSubs(Double_t psi_east, Double_t psi_west);
    void FillPsiRecSubs_wt(Double_t psi_east, Double_t psi_west);
    void FillPsiRecSubs_tpc(Double_t psi_east, Double_t psi_west);
    
    void FillPsiShift(Int_t isub, Double_t psi);
    void FillPsiShiftFull( Double_t psi);
    void FillPsiShiftSubs( Double_t psi0, Double_t psi1);
   
    void FillPsiShift_wt(Int_t isub, Double_t psi);
    void FillPsiShiftFull_wt( Double_t psi);
    void FillPsiShiftSubs_wt( Double_t psi0, Double_t psi1);
    
    void FillPsiShift_tpc(Int_t isub, Double_t psi);
    void FillPsiShiftFull_tpc( Double_t psi);
    void FillPsiShiftSubs_tpc( Double_t psi0, Double_t psi1);

    void FillSubEpQvecRec(Int_t isub, Int_t centnumber, Int_t runindex, Double_t qx, Double_t qy);
    void FillSubEpShiftpar(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    void FillFullEpShiftpar(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    void FillSubEpShiftpar_wt(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    void FillFullEpShiftpar_wt(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    void FillSubEpShiftpar_tpc(Int_t isub, Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    void FillFullEpShiftpar_tpc(Int_t iorder, Int_t centnumber, Int_t runindex, Double_t psi);
    
    void FillPsi1ResolutionEpd(Int_t centnumber, Double_t psi0, Double_t psi1);
    void FillPsi1ResolutionEpdWt(Int_t centnumber, Double_t psi0, Double_t psi1);
    void FillPsi2ResolutionTpc(Int_t centnumber, Double_t psi0, Double_t psi1);
    Double_t GetEtaWeight(Int_t centrality, Double_t eta);
    /*void FillEventEast_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventWest_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j, Int_t k); // i = vertex pos/neg, j = eta_gap, k = ShiftOrder
    void FillEventFull_EP(TVector2 Psi2Vector, TVector2 Psi3Vector, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t k); // i = vertex pos/neg, k = ShiftOrder*/

    // Scalor Product method

  private:

    //raw EventPlane distribution (before recenter and shift)

    // ReCenter Correction | x axis is RunIndex, y axis is Centrality
    //TH1F* href;
    //TH1F* hcent;
  
    TH2D *h2_RingIDvsEta_center;
    TH2D *h2_RingIDvsEta_random;
  
    TH1D *h_phi_center;
    TH1D *h_eta_center;
    TH1D *h_phi_random;
    TH1D *h_eta_random;
  
    TProfile2D *p_sector_eta;
  
    TH1F *h_eta[9];
    TH2D *h2_eta_pt[9];
  
    TProfile2D *p_EPDQx[7];
    TProfile2D *p_EPDQy[7];
    TProfile2D *p_TPCqx_A;
    TProfile2D *p_TPCqy_A;
    TProfile2D *p_TPCqx_B;
    TProfile2D *p_TPCqy_B;
    TString mout_recenter;
    //eta weight
    TH2D *v1EtaWt;
    
    TH1F *h_psi1_epd_ABCD_raw_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_raw_subs;
    TH1F *h_psi1_epd_ABCD_raw_full;
    TH1F *h_psi1_epd_ABCD_recentered_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_recentered_subs;
    TH1F *h_psi1_epd_ABCD_recentered_full;
    TH1F *h_psi1_epd_ABCD_shifted_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_shifted_subs;
    TH1F *h_psi1_epd_ABCD_shifted_full;

    TH1F *hQxSub[_numSubEvents];
    TH1F *hQySub[_numSubEvents];
    TH2F *hQVectorSub[_numSubEvents];
    TH1F *hQxSubRec[_numSubEvents];
    TH1F *hQySubRec[_numSubEvents];
    TH2F *hQVectorSubRec[_numSubEvents];
    TProfile2D *p_mq1x_epd_ABCD_sub[_numSubEvents];
    TProfile2D *p_mq1y_epd_ABCD_sub[_numSubEvents];
    TProfile2D *p_mq1x_epd_ABCD_full;
    TProfile2D *p_mq1y_epd_ABCD_full;
    
    TH1F *h_psi1_epd_ABCD_raw_wt_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_raw_wt_subs;
    TH1F *h_psi1_epd_ABCD_raw_wt_full;
    TH1F *h_psi1_epd_ABCD_recentered_wt_full;
    TH1F *h_psi1_epd_ABCD_recentered_wt_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_recentered_wt_subs;
    TH1F *h_psi1_epd_ABCD_shifted_wt_sub[_numSubEvents];
    TH2F *h_psi1_epd_ABCD_shifted_wt_subs;
    TH1F *h_psi1_epd_ABCD_shifted_wt_full;
    TProfile2D *p_mq1x_epd_ABCD_wt_sub[_numSubEvents];
    TProfile2D *p_mq1y_epd_ABCD_wt_sub[_numSubEvents];
    TProfile2D *p_mq1x_epd_ABCD_wt_full;
    TProfile2D *p_mq1y_epd_ABCD_wt_full;
    
    // TPC event plane
    TH1F *h_psi2_tpc_AB_raw_sub[_numSubEvents];
    TH2F *h_psi2_tpc_AB_raw_subs;
    TH1F *h_psi2_tpc_AB_raw_full;
    TH1F *h_psi2_tpc_AB_recentered_full;
    TH1F *h_psi2_tpc_AB_recentered_sub[_numSubEvents];
    TH2F *h_psi2_tpc_AB_recentered_subs;
    TH1F *h_psi2_tpc_AB_shifted_full;
    TH1F *h_psi2_tpc_AB_shifted_sub[_numSubEvents];
    TH2F *h_psi2_tpc_AB_shifted_subs;
    TProfile2D *p_mq2x_tpc_AB_sub[_numSubEvents];
    TProfile2D *p_mq2y_tpc_AB_sub[_numSubEvents];
    TProfile2D *p_mq2x_tpc_AB_full;
    TProfile2D *p_mq2y_tpc_AB_full;

    // Event Plane method, shift parameter p_TPCshiftpar_Bsin
    TProfile3D *p_shiftpar_sin_epd_ABCD_sub[_numSubEvents];
    TProfile3D *p_shiftpar_cos_epd_ABCD_sub[_numSubEvents];
    TProfile3D *p_shiftpar_sin_epd_ABCD_full;
    TProfile3D *p_shiftpar_cos_epd_ABCD_full;

    TProfile3D *p_shiftpar_sin_epd_ABCD_wt_sub[_numSubEvents];
    TProfile3D *p_shiftpar_cos_epd_ABCD_wt_sub[_numSubEvents];
    TProfile3D *p_shiftpar_sin_epd_ABCD_wt_full;
    TProfile3D *p_shiftpar_cos_epd_ABCD_wt_full;

    TProfile3D *p_shiftpar_sin_tpc_AB_sub[_numSubEvents];
    TProfile3D *p_shiftpar_cos_tpc_AB_sub[_numSubEvents];
    TProfile3D *p_shiftpar_sin_tpc_AB_full;
    TProfile3D *p_shiftpar_cos_tpc_AB_full;
    
    TProfile *p_r1_epd_ABCD_sub0_1, *p_r1_epd_ABCD_wt_sub0_1, *p_r2_tpc_AB_sub0_1;
    TProfile *p_r1_epd_ABCD_sub0_1_sin, *p_r1_epd_ABCD_wt_sub0_1_sin, *p_r2_tpc_AB_sub0_1_sin;
    /*TProfile2D *p_mq2y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq2x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq2y_Full_EP[2];    // 0 = vertex pos/neg

    TProfile2D *p_mq3x_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_East_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3y_West_EP[2][4]; // 0 = vertex pos/neg, 1 = eta gap
    TProfile2D *p_mq3x_Full_EP[2];    // 0 = vertex pos/neg
    TProfile2D *p_mq3y_Full_EP[2];    // 0 = vertex pos/neg*/

    // Shift Correction | x axis is RunIndex, y axis is Centrality
    // Event Plane method
    /*TProfile2D *p_mcos2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin2_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin2_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder

    TProfile2D *p_mcos3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_East_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_msin3_West_EP[2][4][5]; // 0 = vertex pos/neg, 1 = eta gap, 2 = ShiftOrder
    TProfile2D *p_mcos3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder
    TProfile2D *p_msin3_Full_EP[2][5];    // 0 = vertex pos/neg, 1 = ShiftOrder*/


    static TString mVStr[2];

    ClassDef(EpProManager,1)
};

#endif
