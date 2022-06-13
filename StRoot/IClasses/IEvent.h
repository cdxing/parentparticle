// author: Erik Loyd
// Based off of code by Mike Lisa

#ifndef IEVENT_H
#define IEVENT_H

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <vector>
#include "IEventPlane.h"


/*
  Store simple event summary information
*/

class IEvent : public TObject {

 protected:

/*
  Int_t mRunNumber;
  Int_t mChronologicalRunId;
  Int_t mEventID;  // within run
  Float_t mBfield;
  TVector3 mPrimaryVertex;
  Int_t mRefMult;
  Int_t mCentrality;
  Float_t mVpdVz;
  std::vector<unsigned int> mTriggers;
  Short_t             EPDid[744];
  UShort_t            nEPDhits;

  UShort_t            tree_tracknumber; //What is this?
 
*/


  Float_t EPDnMip[2][12][31]; //0 = E, 1 = W
 
  Float_t mqCenter[2];  // This is the Q recentering offset for all events in this data set
  std::vector<IEventPlane> mEPParticles;
  
  void CalcQVector(int, std::vector<float>&, std::vector<float>&);
  std::vector<float> CalcPsi(int);
 
 public:

  IEvent();
  ~IEvent();

  void ClearEvent();
  void Init();

  //Vector with event psi as well as auto-correlation-removed psis
  std::vector<float> GetPsi(int);
  std::vector<float> GetPsi(int, char, float, float);  
  
  //Get the Q components for the event
  Float_t GetQx(int);
  Float_t GetQy(int);
  
  //Get the Psi for the event (or for a subevent)
  Float_t GetEventPsi(int harmonic){return GetPsi(harmonic)[0];}
  Float_t GetEventPsi(int harmonic, char option, float param1, float param2 = 0.0){return GetPsi(harmonic, option, param1, param2)[0];}
  
  //Get the autocorrelation-corrected Psi for the each particle in an event (or subevent)
  std::vector<float> GetAutoCorrelationPsi(int);
  std::vector<float> GetAutoCorrelationPsi(int, char, float, float);

  std::vector<IEventPlane> GetEPParticles(){return mEPParticles;}
  
  void AddEPParticle(IEventPlane n){mEPParticles.push_back(n);} //Add a particle to the event
  void AddEPParticle(std::vector<IEventPlane> n){mEPParticles = n;} //Add a vector of particles to the event
  
  void RemoveEPParticle(int n){mEPParticles.erase(mEPParticles.begin()+n);} //Remove a particle at index n


  void SetEPParticles( const std::vector<IEventPlane> & n){mEPParticles = n;} //Set the particles for the event
  void SetQCenter(float, float); //Set the center for the Q vector - for recentering
  
  std::vector<IEventPlane> EPDVector(TVector3, float, float, float);
  
  void SetEPDnMip(int ew, int pp, int tt, float nMip){EPDnMip[ew][pp - 1][tt - 1] = nMip;}
  
  void SetEPDnMip(int tileID, float nMip);
  
  Float_t GetEPDnMip(int ew, int pp, int tt){return EPDnMip[ew][pp - 1][tt - 1];}
  Float_t GetEPDnMip(int tileID);

  IEvent GetSubEvent(char, float, float);
  
  void AddEPDtoTracks(TVector3, float);

  ClassDef(IEvent,1)  // my event
};

#endif