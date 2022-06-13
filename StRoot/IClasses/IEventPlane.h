// Tiny DST class for calculating event planes
#ifndef IEVENTPLANE_H
#define IEVENTPLANE_H

#include "TObject.h"
#include "TVector3.h"
#include "TMath.h"

class IEventPlane : public TObject {

	private:

		Float_t mPhi;
		Float_t mWeight;

		Float_t mEta;
		TVector3 mMomentum;
		
		Int_t mPID;
		Int_t mCharge;
		Float_t mToFBeta;
		Int_t nHitsFit;
		Int_t nHitsPoss;
		Float_t dEdx;
		
		Int_t mtileID;
		Float_t nMIP;
		
		void Construct();


	public:
		IEventPlane();
		IEventPlane(float, float);
		
		~IEventPlane(){} //This is needed apparently

		void SetPhi(){mPhi = mMomentum.Phi();}
		void SetPhi(float in_phi){mPhi = in_phi;}
		
		void SetEta(){mEta = mMomentum.Eta();}
		void SetEta(float in_eta){mEta = in_eta;}		
		
		void SetWeight(float in_weight){mWeight = in_weight;}
		
		void SetMomentum(TVector3 in_momentum){mMomentum = in_momentum;}
		void SetMomentumX(float in_X){mMomentum.SetX((double)in_X);}
		void SetMomentumY(float in_Y){mMomentum.SetY((double)in_Y);}
		void SetMomentumZ(float in_Z){mMomentum.SetZ((double)in_Z);}
		void SetMomentumXYZ(float in_X, float in_Y, float in_Z){mMomentum.SetXYZ(in_X,in_Y,in_Z);}
		
		void SetCharge(int in_Charge){mCharge = in_Charge;}
		void SetHitsFit(float in_Fit){nHitsFit = in_Fit;}
		void SetHitsPoss(float in_Poss){nHitsFit = in_Poss;}
		void SetParticleID(int in_PID){mPID = in_PID;}
		void SetToFBeta(float in_ToF){mToFBeta = in_ToF;}
		void SetTileID(int in_TileID){mtileID = in_TileID;}
		void SetnMIP(float in_mip){nMIP = in_mip;}
		void SetdEdx(float in_dedx){dEdx = in_dedx;}
	

		Float_t GetPhi(){return mPhi;}
		Float_t GetWeight(){return mWeight;}
		
		Int_t GetCharge(){return mCharge;}		
		Float_t GetEta(){return mEta;}		
		Int_t GetHitsFit(){return nHitsFit;}
		Int_t GetHitsPoss(){return nHitsPoss;}		
		TVector3 GetMomentum(){return mMomentum;}
		Int_t GetParticleID(){return mPID;}				
		Float_t GetToFBeta(){return mToFBeta;}
		Float_t GetdEdx(){return dEdx;}
		
		
		Int_t GetEPDtileID(){return mtileID;}
		Int_t GetEPDew();
		Int_t GetEPDtt(){return TMath::Abs(mtileID) % 100;}
		Int_t GetEPDpp(){return (TMath::Abs(mtileID) - GetEPDtt()) / 100;}
		Int_t GetRingNumber(){return TMath::Ceil((GetEPDtt() + 1) / 2);}
		
		Float_t	GetnMIP(){return nMIP;}
		
		
		Float_t QxTerm(int);
		Float_t QyTerm(int);

		ClassDef(IEventPlane,1)  // my EventPlane
};
    
#endif