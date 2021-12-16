
clc,clear all,close all

load ('D:\Tasks\BasicSWEController\DTU10MW_FBSWE2_DLC1d2\DTU10MW_FBSWE2_DLC1d2_statistics.mat');
cases ={ 'V4_TSR72';
         'V4_TSR74';
         'V4_TSR76';
         'V4_TSR78';
         'V6_TSR72';
         'V6_TSR74';
         'V6_TSR76';
         'V6_TSR78';
         'V8_TSR72';
         'V8_TSR74';
         'V8_TSR76';
         'V8_TSR78';
           };
       
Pelmean =  Statistics.MEAN(:,8);
MytDel  =  Statistics.DEL(:,6);
LSSDel  =  Statistics.DEL(:,9);
MyBlDel =  Statistics.DEL(:,10);   
   
MytDel4    =  Statistics.DEL(1:4,6);  
MytDel4rel =  100*(MytDel4./ MytDel4(1))-100;
MytDel6    =  Statistics.DEL(5:8,6);
MytDel6rel =  100*(MytDel6./ MytDel6(1))-100;
MytDel8    =  Statistics.DEL(9:12,6);
MytDel8rel = 100*( MytDel8./ MytDel8(1))-100;
MytDel10    =  Statistics.DEL(13:16,6);
MytDel10rel = 100*( MytDel10./ MytDel10(1))-100;

Pelmean4  = Pelmean(1:4);
Pelmean4rel =  100*(Pelmean4./ Pelmean4(1))-100;
Pelmean6  = Pelmean(5:8);
Pelmean6rel =  100*((Pelmean6./ Pelmean6(1)))-100;
Pelmean8  = Pelmean(9:12);
Pelmean8rel =  100*((Pelmean8./ Pelmean8(1)))-100;
Pelmean10    =  Pelmean(13:16);
Pelmean10rel = 100*( Pelmean10./ Pelmean10(1))-100;

MyBlDel4    =  MyBlDel(1:4);  
MyBlDel4rel =  100*(MyBlDel4./ MyBlDel4(1))-100;
MyBlDel6    =  MyBlDel(5:8);
MyBlDel6rel =  100*(MyBlDel6./ MyBlDel6(1))-100;
MyBlDel8    =  MyBlDel(9:12);
MyBlDel8rel = 100*( MyBlDel8./ MyBlDel8(1))-100;
MyBlDel10    =  Pelmean(13:16);
MyBlDel10rel = 100*( MyBlDel10./ MyBlDel10(1))-100;

LSSDel4    =  LSSDel(1:4);  
LSSDel4rel =  100*(LSSDel4./ LSSDel4(1))-100;
LSSDel6    =  LSSDel(5:8);
LSSDel6rel =  100*(LSSDel6./ LSSDel6(1))-100;
LSSDel8    =  LSSDel(9:12);
LSSDel8rel = 100*( LSSDel8./ LSSDel8(1))-100;
LSSDel10    =  Pelmean(13:16);
LSSDel10rel = 100*( LSSDel10./ LSSDel10(1))-100;









