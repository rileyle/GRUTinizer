#ifndef TRAWGSBANKS_H__
#define TRAWGSBANKS_H__


#include <Rtypes.h>



struct GS_TapeHeader    {
  UShort_t RecordType;      /* tape header type= 1             */
  UShort_t RecordLength;    /* number of bytes in this record  */
  UShort_t RecordVer;       /* record version or subtype       */
  UInt_t   ByteOrder;       /* to determine byte ordering      */
  char     exp_title1[40];  /* first part of exp title         */
  char     exp_title2[40];  /* 2nd  part of exp title          */
  char     time[9];         /* hh:mm:ss         */
  char     date[9];         /* yy/mm/dd         */
  UShort_t TapeNum;         /* tape number in this experiment  */
  UShort_t TapeUnit;        /* specifies which unit wrote this tape */

  void swap() {
    RecordType   = __builtin_bswap16(RecordType);        
    RecordLength = __builtin_bswap16(RecordLength);    
    RecordVer    = __builtin_bswap16(RecordVer);       
    ByteOrder    = __builtin_bswap32(ByteOrder);       
    TapeNum      = __builtin_bswap16(TapeNum);       
    TapeUnit     = __builtin_bswap16(TapeUnit);       
  }

  void print() const {
    printf("RecordType:  %12i\n",RecordType);
    printf("RecordLength:%12i\n",RecordLength);
    printf("RecordVer:   %12i\n",RecordVer);
    printf("ByteOrder:    0x%08x\n",ByteOrder);
    printf("ExpTitle1:   %12s\n",exp_title1);
    printf("ExpTitle2:   %12s\n",exp_title2);
    printf("Time:        %12s\n",time);
    printf("Date:        %12s\n",date);
    //return;
    printf("TapeNum:           0x%04x\n",TapeNum);
    printf("Tapeunit:          0x%04x\n",TapeUnit);
  }
}__attribute__((__packed__)); 


struct GS_FileHeader {
  UShort_t RecordType;     /* file header type= 2             */
  UShort_t RecordLength;   /* number of bytes in this record  */
  UShort_t RecordVer;      /* record version or subtype       */
  UShort_t RunNumber;      /* run number in this experiment   */
  UShort_t FileNumber;     /* file number in this experiment  */
  char    run_title1[40];  /* first part of run title         */
  char    run_title2[40];  /* 2nd  part of run title          */

  void swap() {
    RecordType   = __builtin_bswap16( RecordType);  
    RecordLength = __builtin_bswap16(RecordLength);
    RecordVer    = __builtin_bswap16(RecordVer);   
    RunNumber    = __builtin_bswap16(RunNumber);   
    FileNumber   = __builtin_bswap16(FileNumber);   
  };

  void print() {
    printf("RecordType:    % 8i\n",RecordType);    
    printf("RecordLength:  % 8i\n",RecordLength);
    printf("RecordVer:     % 8i\n",RecordVer);  
    printf("RunNumber:     % 8i\n",RunNumber);  
    printf("FileNumber:    % 8i\n",FileNumber);      
    printf("run_title1:    %8s\n",run_title1);
    printf("run_title2:    %8s\n",run_title2);
  };
}__attribute__((__packed__));   


struct GS_EventBuffer  {
  //UShort_t RecordType;        /* 0 event data type= 3              */
  //UShort_t RecordLength;      /* 1 number of bytes in this record  */
  UShort_t RecordVer;         /* 2 record version or subtype       */
  UShort_t HeaderBytes;       /* 3 number of bytes in header       */
  UShort_t EffNumber;         /* 4 eff processor number            */
  UShort_t StreamID;          /* 5 event stream ID                 */
  UShort_t EffSequence;       /* 6 eff sequence number             */
  UShort_t ModeFlags;         /* 7 event format flags              */
  UShort_t DataLength;        /* 8 number of i*2 data words        */
  UShort_t ChecksumType;      /* 9 type of checksum                */
  UShort_t Checksum;          /* 10 checksum value                 */
  //u_short EventData[EB_SIZE]; /* event data area                   */

  void swap() {
    //RecordType   = __builtin_bswap16(RecordType);  
    //RecordLength = __builtin_bswap16(RecordLength);
    RecordVer    = __builtin_bswap16(RecordVer);  
    HeaderBytes  = __builtin_bswap16(HeaderBytes);  
    EffNumber    = __builtin_bswap16(EffNumber);  
    StreamID     = __builtin_bswap16(StreamID);  
    EffSequence  = __builtin_bswap16(EffSequence);  
    ModeFlags    = __builtin_bswap16(ModeFlags);  
    DataLength   = __builtin_bswap16(DataLength);  
    ChecksumType = __builtin_bswap16(ChecksumType);
    Checksum     = __builtin_bswap16(Checksum);  
  };

  void print() {
    //printf("RecordType:     % 8i\n",RecordType);  
    //printf("RecordLength:   % 8i\n",RecordLength);
    printf("RecordVer:      % 8i\n",RecordVer);  
    printf("HeaderBytes:    % 8i\n",HeaderBytes);  
    printf("EffNumber:      % 8i\n",EffNumber);  
    printf("StreamID:       % 8i\n",StreamID);  
    printf("EffSequence:    % 8i\n",EffSequence);  
    printf("ModeFlags:        0x%04x\n",ModeFlags);  
    printf("DataLength:     % 8i\n",DataLength);  
    printf("ChecksumType:   % 8i\n",ChecksumType);
    printf("Checksum:         0x%04x\n",Checksum);  
  };

}__attribute__((__packed__));   


struct GS_EventHeader {
  UShort_t word0;  // 0x80nn    total #words in event         
  UShort_t word1;  // 0xmmnn    nn clean ge counts; mm is the MTM bit pattern at MAIN
  UShort_t word2;  // 0xmmnn    dirty ge and bgo only counts, respectively
  UShort_t word3;  // 0xffff    trigger time, ttH (16 bits)
  UShort_t word4;  // 0x7fff    trigger time, ttM (15 bits)
  UShort_t word5;  // 0xffff    trigger time, ttL (16 bits)
  UShort_t word6;  // 0x0nnn    tac1(usec tac; next usec tic vs. master trigger)
  UShort_t word7;  // 0x0nnn    tac2 (usually rf tac: next rf vs. pre-trigger)
  UShort_t word8;  // 0xnnnn    sum of ge_low energies (within ge time window)
  UShort_t word9;  // 0xnnnn    sum of BGO energies    (within BGO time window)

  void swap() {
    word0 = __builtin_bswap16(word0);  
    word1 = __builtin_bswap16(word1);
    word2 = __builtin_bswap16(word2);
    word3 = __builtin_bswap16(word3);
    word4 = __builtin_bswap16(word4);
    word5 = __builtin_bswap16(word5);
    word6 = __builtin_bswap16(word6);
    word7 = __builtin_bswap16(word7);
    word8 = __builtin_bswap16(word8);
    word9 = __builtin_bswap16(word9);
  }

  UShort_t total_words() { return word0&0x00ff; }
  UShort_t clean_hpge()  { return word1&0x00ff; }
  UShort_t mtm_bitp()    { return (word1&0xff00)>>8; }
  UShort_t dirty_hpge()  { return (word2&0xff00)>>8; }
  UShort_t bgo_only()    { return (word2&0x00ff); }
  uint64_t trig_time()   { return (((uint64_t)word3)<<31) + (((uint64_t)word4)<<16) + (((uint64_t)word5)<<0); }
  UShort_t tac1()        { return word6&0xfff; }
  UShort_t tac2()        { return word7&0xfff; }
  UShort_t hpge_sum()    { return word8;   }
  UShort_t bgo_sum()     { return word9;   }

  void print() {
    printf("total_words:  % 6i\n", total_words());
    printf("clean_hpge:   % 6i\n", clean_hpge() );
    printf("mtm_bitp:     % 6i\n", mtm_bitp()   );
    printf("dirty_hpge:   % 6i\n", dirty_hpge() );
    printf("bgo_only:     % 6i\n", bgo_only()   );
    printf("trig_time:    %lu\n", trig_time()  );
    printf("tac1:         % 6i\n", tac1()       );
    printf("tac2:         % 6i\n", tac2()       );
    printf("hpge_sum:     % 6i\n", hpge_sum()   );
    printf("bgo_sum:      % 6i\n", bgo_sum()    );
  }

}__attribute__((__packed__));





struct GS_HpgeClean {
  UShort_t hpid;
  UShort_t ge_high;
  UShort_t ge_side;
  UShort_t ge_time;
  UShort_t ge_trap;
  UShort_t ge_low;

  void swap() {
    hpid    = __builtin_bswap16(hpid);
    ge_high = __builtin_bswap16(ge_high);
    ge_side = __builtin_bswap16(ge_side);
    ge_time = __builtin_bswap16(ge_time);
    ge_trap = __builtin_bswap16(ge_trap);
    ge_low  = __builtin_bswap16(ge_low);
  }

  UShort_t bgo_hitp()   { return (hpid&0xfe00) >> 9; }
  bool     ge_hit_bit() { return (hpid&0x0100) >> 8; }
  UShort_t id()         { return (hpid&0x00ff) >> 0; }
  UShort_t charge()     { return (ge_high&0x3fff);   }
  bool     over_range() { return (ge_high&0x4000)>>14;}
  bool     pileup()     { return (ge_high&0x8000)>>15;}
  UShort_t side_chg()   { return (ge_side&0x0fff);   }
  UShort_t time()       { return (ge_time&0x1fff);   }

  void print() { printf("hpge[%03i] = % 8i  | % 8i\n",id(),charge(),time()); }



}__attribute__((__packed__));

struct GS_HpgeDirty {
  UShort_t hpid;
  UShort_t ge_high;
  UShort_t ge_side;
  UShort_t ge_time;
  UShort_t ge_trap;
  UShort_t ge_low;
  UShort_t bgo_time;
  UShort_t bgo_low;


  void swap() {
    hpid    = __builtin_bswap16(hpid);
    ge_high = __builtin_bswap16(ge_high);
    ge_side = __builtin_bswap16(ge_side);
    ge_time = __builtin_bswap16(ge_time);
    ge_trap = __builtin_bswap16(ge_trap);
    ge_low  = __builtin_bswap16(ge_low);
    bgo_time = __builtin_bswap16(bgo_time);
    bgo_low = __builtin_bswap16(bgo_low);    
  }

  UShort_t bgo_hitp()   { return (hpid&0xfe00) >> 9; }
  bool     ge_hit_bit() { return (hpid&0x0100) >> 8; }
  UShort_t id()         { return (hpid&0x00ff) >> 0; }
  UShort_t charge()     { return (ge_high&0x3fff);   }
  bool     over_range() { return (ge_high&0x4000)>>14;}
  bool     pileup()     { return (ge_high&0x8000)>>15;}
  UShort_t side_chg()   { return (ge_side&0x0fff);   }
  UShort_t time()       { return (ge_time&0x1fff);   }
  UShort_t btime()      { return (bgo_time&0x0fff);  }
  UShort_t bgo_chg()    { return (bgo_low&0x0fff);   }

  void print() { printf("hpge[%03i] = % 8i  | % 8i  bgo: % 8i | % 8i\n",id(),charge(),time(),
      bgo_chg(),btime()); }

}__attribute__((__packed__));

struct GS_BgoClean {
  UShort_t hpid;
  UShort_t bgo_time;
  UShort_t bgo_low;

  void swap() {
    hpid    = __builtin_bswap16(hpid);
    bgo_time = __builtin_bswap16(bgo_time);
    bgo_low = __builtin_bswap16(bgo_low);    
  }

  UShort_t bgo_hitp()   { return (hpid&0xfe00) >> 9; }
  bool     ge_hit_bit() { return (hpid&0x0100) >> 8; }
  UShort_t id()         { return (hpid&0x00ff) >> 0; }
  UShort_t charge()     { return (bgo_low&0x0fff);   }
  UShort_t time()       { return (bgo_time&0x0fff);   }

  void print() { printf("bgo[%03i] = % 8i  | % 8i\n",id(),charge(),time()); }

}__attribute__((__packed__));





#endif


