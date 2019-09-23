#ifndef _TOBJHIT_H_
#define _TOBJHIT_H_

#include <TVector3.h>
#include <TDetectorHit.h>
#include <GH1D.h>

class TOBJHit : public TDetectorHit {
  public:
    TOBJHit() : htrace(nullptr),pileup(false) { }
    ~TOBJHit() { }

    virtual void Copy(TObject& obj) const        { TDetectorHit::Copy(obj); }
    virtual void Clear(Option_t *opt = "");//    { TDetectorHit::Clear(opt); }
    virtual void Print(Option_t *opt = "") const { TDetectorHit::Print(opt); }
    //virtual void Draw(Option_t* opt = "")        { TDetectorHit::Draw(opt); }

    std::vector<unsigned short>* GetTrace()   { return &fTrace; } 
    void SetTrace(unsigned int trace_length, const unsigned short* trace);
    //void DrawTrace() { }
    int Size() { return fTrace.size(); }

    int GetCrate()   const { return 0; }
    int GetSlot()    const { return ((Address()&0x000000f0)>>4); }
    int GetChannel() const { return ((GetSlot()-2)*16) + (Address()&0x0000000f); }

    TVector3 GetPosition() { return TVector3(0,0,1); } 
   
    int GetNumber() const  { return GetSlot()*16+GetChannel(); }

    void SetExternalTimestamp(long ets) { fExternalTimestamp = ets; }
    long GetExternalTimestamp() const { return fExternalTimestamp; }
    void DrawTrace(Option_t *option);
    void TrigFilter(Double_t tpeak,Double_t tgap); // in unit of ticks
    void EnergyFilter(Double_t tpeak,Double_t tgap,Double_t tau); // in unit of ticks
    void SetPileup(bool flag){pileup = flag; }
    bool GetPileip() const { return pileup;}
    
 

  private:
    std::vector<unsigned short> fTrace;
    std::vector<unsigned short> fETrace;
    std::vector<unsigned short> fTTrace;
    GH1D *htrace;
    Bool_t  pileup;

    long fExternalTimestamp;

  ClassDef(TOBJHit,1)
};

#endif
