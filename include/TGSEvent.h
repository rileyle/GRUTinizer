#ifndef _TGSEVENT_H_
#define _TGSEVENT_H_

#include "TRawEvent.h"

class TMode3;

class TGSEvent : public TRawEvent {

public:
#include "TRawGSBanks.h"

  public:
    TGSEvent();
    TGSEvent(const TRawEvent&);
    virtual ~TGSEvent();

    enum EEventType {kTapeHeader=1,kFileHeader=2,kDataRecord=3 };


    long GetTimestamp() const;

    const char* GetPayload() const;
    TSmartBuffer GetPayloadBuffer() const;

    virtual void Clear(Option_t *opt ="");
    virtual void Print(Option_t *opt ="") const;

  ClassDef(TGSEvent, 0);
};




#endif /* _TGSEVENT_H_ */
