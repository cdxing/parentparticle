#include "TMath.h"
#include "IEventPlane.h"

ClassImp(IEventPlane)

//______________
IEventPlane::IEventPlane(){
	mPhi = 0;
	SetWeight(0);
	Construct();
}

//______________
IEventPlane::IEventPlane(float in_phi, float in_weight){
	mPhi = in_phi;
	SetWeight(in_weight);
	Construct();
}

void IEventPlane::Construct(){
	mEta = 0.0;
	mMomentum.SetXYZ(0,0,0);
	mPID = 0;
	mCharge = 0;
	mToFBeta = 0.0;
	nHitsFit = 0;
	nHitsPoss = 0;
	mtileID = 0;
	nMIP = 0.0;
}

//______________
Float_t IEventPlane::QxTerm(int harmonic){
	Float_t Cosine = cos(mPhi*(float)harmonic);
	return mWeight * Cosine;
}

//______________
Float_t IEventPlane::QyTerm(int harmonic){
	Float_t Sine = sin(mPhi*(float)harmonic);
	return mWeight * Sine;
}

//0 = east, 1 = west
Int_t IEventPlane::GetEPDew(){
	if (mtileID > 0){return 1;}
	return 0;
}