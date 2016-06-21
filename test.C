
{
 TRawFileIn infile("/data0/bender/gs/34Si_A.r001.dat.000.gz",kFileType::GAMMASPHERE_DAT);
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
       //event.Print();
       break;
   }
   //if(counter>10) break;
 }

 printf("--------------------\n");
 for(int x=1;x<4;x++) {
    printf("\ttype[%i]  = %i\n",x,types[x]);
 }
 printf("--------------------\n");


}



