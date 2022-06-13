#ifndef Shift_h
#define Shift_h

#include "StMaker.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/EpProManager/EpProManager.h" // Load STARLibrary header files
#include "TString.h"
#include "TFile.h"
// Configuration file reader
#include "../ConfigReader/ConfigReader.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class CutManager;
class HistManager;
class StRefMultCorr;
class StPicoBTofPidTraits;
class EpProManager;
class FlowCorrection;
class FlowV0;

// run QA
//class TProfile;
//class TH1F;
//class TH2F;

class ParentParticle : public StMaker {
    public:
        ParentParticle(const char *name, StPicoDstMaker *picoMaker, char* jobid, std::string configFileName);
        virtual ~ParentParticle();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        //int Centrality(int gRefMult);
        int GetRunIndex(int runId);

        //bool isGoodTrigger(StPicoEvent const*) const;

    private:
	ConfigReader configs;
        static StRefMultCorr *mRefMultCorr;
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *mPicoEvent;
        CutManager  *mCutManager;
        FlowCorrection *mFlowCorrection;
        FlowV0 *mFlowV0;
        //HistManager *mHistManager;

	// EPD 
        StEpdGeom *mEpdGeom;
        StEpdEpInfo *mEpdEpInfo;
        EpProManager *mEpProManager;

	Int_t mEnergy;

	// get the recenterpar
        TProfile2D *getep_sub_recen[2][_numSubEvents];
        TProfile2D *getep_full_recen[2];

        TProfile2D *getep_sub_wt_recen[2][_numSubEvents];
        TProfile2D *getep_full_wt_recen[2];

        TProfile2D *getep_sub_tpc_recen[2][_numSubEvents];
        TProfile2D *getep_full_tpc_recen[2];
        
	//  get shift par
        TProfile3D *p_sub_ep_shiftpar_sin[_numSubEvents];
        TProfile3D *p_sub_ep_shiftpar_cos[_numSubEvents];
        TProfile3D *p_full_ep_shiftpar_sin;
        TProfile3D *p_full_ep_shiftpar_cos;
        
        TProfile3D *p_sub_ep_wt_shiftpar_sin[_numSubEvents];
        TProfile3D *p_sub_ep_wt_shiftpar_cos[_numSubEvents];
        TProfile3D *p_full_ep_wt_shiftpar_sin;
        TProfile3D *p_full_ep_wt_shiftpar_cos;

        TProfile3D *p_sub_ep_tpc_shiftpar_sin[_numSubEvents];
        TProfile3D *p_sub_ep_tpc_shiftpar_cos[_numSubEvents];
        TProfile3D *p_full_ep_tpc_shiftpar_sin;
        TProfile3D *p_full_ep_tpc_shiftpar_cos;
	TString mOutPut_Phi;

        TFile *mFile_Phi;


        ClassDef(ParentParticle, 1)
};

#endif
