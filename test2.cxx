//g++ test2.cxx `grutinizer-config --cflags --all-libs --root` -lTGSFormat -O3


#include <cstdio>

#include <TRawSource.h>
#include <TGSEvent.h>


#include <TFile.h>
#include <TH2.h>

TH2F *hist = new TH2F("hist","hist",120,0,120,8000,0,16000);

int HandleEvent(char *data) {
  int ptr = 0;
  TGSEvent::GS_EventHeader *header = (TGSEvent::GS_EventHeader*)data;
  if(header->total_words()>0xfff || header->total_words()==0)
    return -1;
  header->swap();
  //header->print();
  ptr += sizeof(TGSEvent::GS_EventHeader);
  for(int x=0;x<header->clean_hpge();x++) {  
    TGSEvent::GS_HpgeClean *hpge = (TGSEvent::GS_HpgeClean*)(data+ptr);
    hpge->swap(); 
    //printf("\t"); hpge->print();
    ptr += sizeof(TGSEvent::GS_HpgeClean);
    hist->Fill(hpge->id(),hpge->charge());
  }
  for(int x=0;x<header->dirty_hpge();x++) {  
    TGSEvent::GS_HpgeDirty *hpge = (TGSEvent::GS_HpgeDirty*)(data+ptr);
    hpge->swap(); 
    //printf("\t"); hpge->print();
    ptr += sizeof(TGSEvent::GS_HpgeDirty);
    hist->Fill(hpge->id(),hpge->charge());
  }
  ptr += header->bgo_only()*sizeof(TGSEvent::GS_BgoClean);
 
  //printf("  ptr = %i\n  ",ptr); 
  //printf("--------------------\n");
  //return sizeof(TGSEvent::GS_EventHeader) + header->total_words()*2;
  ptr += 2; // 0xffff  event separator...
  return ptr;
}

bool HandleDataRecord(TGSEvent *event) {
  int ptr=0;;
  TGSEvent::GS_EventBuffer *dataheader;
  dataheader = (TGSEvent::GS_EventBuffer*)event->GetPayload();
  dataheader->swap();
  //dataheader->print();
  ptr+=sizeof(TGSEvent::GS_EventBuffer);
  while(ptr<dataheader->DataLength) {
    ptr += HandleEvent((char*)(event->GetPayload()+ptr));
  }
  //printf("|||||||||||||||||||||||||\n");
  //printf("|||||||||||||||||||||||||\n");
}

int main(int argc, char **argv) {
 TRawFileIn infile("/data0/bender/gs/Eu152.r027.dat.000.gz",kFileType::GAMMASPHERE_DAT);
 TGSEvent event;
 

 int counter=0;
 int types[9] = {0};

 while(infile.Read(&event)>0) { 
   counter++; 
   switch(event.GetEventType()) {
     case TGSEvent::kTapeHeader:
       types[1]++;
       event.Print();
       break;
     case TGSEvent::kFileHeader:
       types[2]++;
       break;
     case TGSEvent::kDataRecord:
       types[3]++;
       HandleDataRecord(&event);
       //event.Print();
       //exit(0);
       break;
   }
   //if(counter>10) break;
   if((counter%1000)==0) {
     printf("  past record %i           \r",counter);
     fflush(stdout);
   }
 }
 printf("  finished processesing %i records.\n",counter);  

 printf("--------------------\n");
 for(int x=1;x<4;x++) {
    printf("\ttype[%i]  = %i\n",x,types[x]);
 }
 printf("--------------------\n");

 TFile f("test.root","recreate");
 hist->Write();
 f.Close();
 return 0;
}



