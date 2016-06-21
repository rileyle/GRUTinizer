#include "TGSEvent.h"

#include "TString.h"


ClassImp(TGSEvent)

TGSEvent::TGSEvent() { }

TGSEvent::TGSEvent(const TRawEvent &raw) {
  raw.Copy(*this);
}

TGSEvent::~TGSEvent() { }

Long_t TGSEvent::GetTimestamp() const {
  return -1;
  //return *((Long_t*)(TRawEvent::GetBody() + 0));
}

const char* TGSEvent::GetPayload() const {
  //return fBody.GetData() + sizeof(Int_t);
  return fBody.GetData();
  //return fBody.GetData()-sizeof(Int_t);
}

TSmartBuffer TGSEvent::GetPayloadBuffer() const {
  if(fTimestamp != -1){
    return fBody;
  } else {
    return fBody.BufferSubset(sizeof(Int_t));
  }
}

void TGSEvent::Clear(Option_t *opt) {
  TRawEvent::Clear(opt);
}

void TGSEvent::Print(Option_t *opt) const {
  TString options(opt);
  std::cout << BLUE << "Type:     \t" << DYELLOW << GetEventType() << RESET_COLOR << std::endl;
  std::cout << BLUE << "Size:     \t" << DYELLOW << GetBodySize()  << RESET_COLOR << std::endl;
  //std::cout << BLUE << "Timestamp:\t" << DYELLOW << GetTimestamp() << BLUE << "  tens of ns" << RESET_COLOR << std::endl;
  if(options.Contains("all")) {
    TString pass_opt("bodyonly");
    pass_opt.Append(options);
    TRawEvent::Print(pass_opt.Data());
  }
}

//ClassImp(TGSMode3Event)
//
//void TGSMode3Event::BuildFragments(){
//  TSmartBuffer buf = fEvent.GetPayloadBuffer();
//  TGSEvent event(fEvent);
// 
//  while(buf.GetSize()){
//    // Read the header and body
//    TRawEvent::GEBMode3Head* header = (TRawEvent::GEBMode3Head*)buf.GetData();
//    TRawEvent::SwapMode3Head(*header);
// 
//    TRawEvent::GEBMode3Data* data = (TRawEvent::GEBMode3Data*)(buf.GetData()+sizeof(TRawEvent::GEBMode3Head));
//    TRawEvent::SwapMode3Data(*data);
// 
//    //header.GetLength() is number of 32-bit values,
//    //   and does not include the 0xaaaaaaaa separator.
//    size_t body_size = header->GetLength()*4 + 4;
// 
//    // Transfer the timestamp and body
//    event.SetFragmentTimestamp(data->GetLed());
//    event.SetData(buf.BufferSubset(0,body_size));
// 
//    // Push a copy into the list
//    fragments.push_back(event);
//    buf.Advance(body_size);
//  }
// 
//  // std::cout << "---------------------------------------------------" << std::endl;
//  // fEvent.Print("all");
//  // std::cout << "NumFragments: " << fragments.size() << std::endl;
//  // for(auto& frag : fragments){
//  //   frag.Print("all");
//  // }
// 
//  // std::cout << "---------------------------------------------------" << std::endl;
//}
