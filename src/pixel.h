#ifndef _PIXEL_H__
#define _PIXEL_H__

class CPixel
{
public :
   int time;
   int cc;
   float flux;
   float weight;
   
   double exact_uxtime;
   double exact_freq_mhz;
   
   CPixel( int _time=0, int _cc=0, float _flux=0.00 ): 
   time(_time), cc(_cc), flux(_flux), exact_uxtime(0), exact_freq_mhz(0), weight(1.00)
   {       
   }
};


#endif
